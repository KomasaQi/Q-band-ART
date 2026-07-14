%% *****************MPC控制器介绍*********************************
% 模型：tractor-side 2DOF single-track dynamics + Q-band ART
% 在线预测只使用牵引车 vy、yaw rate、vx、heading、Y、X。
% 液体状态、挂车侧状态和真实轮载/LTR 不作为 MPC 状态。
% 观测量：observe=[Y1,X1,psi,vy,r,LTR_est]';6×1
% 控制器：轨迹跟踪 + 横向加速度危险频带抑制

%% 主函数 根据标志位flag来切换操作具体方程
% MPC_TrajTrackT1_TPSP_LTR_Ctrllor
function [sys,x0,str,ts] = Pin5_MPC_Qband_TrajTrack_Ctrllor(t,x,u,flag)
% 本测试脚本用于尝试纠正MPC车辆预测模型中的潜在错误
    %   限定于车辆动力学模型，控制量为前轮偏角
    %   [sys,x0,str,ts] = MY_MPCController3(t,x,u,flag)
    switch flag
     case 0
      [sys,x0,str,ts] = mdlInitializeSizes; % Initialization 
     case 2
      sys = mdlUpdates(t,x,u); % Update discrete states
     case 3
      sys = mdlOutputs(t,x,u); % Calculate outputs
     case {1,4,9} % Unused flags
      sys = []; 
     otherwise
      error(['unhandled flag = ',num2str(flag)]); % Error handling
    end
    % End of dsfunc.
end
%% 模块初始化
%==============================================================
% Initialization
%==============================================================
function [sys,x0,str,ts] = mdlInitializeSizes

    % Call simsizes for a sizes structure, fill it in, and convert it 
    % to a sizes array.
    
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 1; % this parameter doesn't matter
    sizes.NumOutputs     = 4; % [delta, M1zp, M1z_, M2zp, M2z_, vdes,]'+
%                               [LTR, epsilon, iter_time];
    sizes.NumInputs      = 16;% TruckSim 仍输入 16 个量；控制器只使用牵引车侧
%                               vy1, yaw rate, vx1, heading, Y1, X1。
%                               其余液体/挂车相关输入仅保留接口兼容，不进入 MPC。
    sizes.DirFeedthrough = 1; % Matrix D is non-empty.
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes); 
    x0 =1e-4;   
    global U; % store current ctrl vector: delta_m
    U=0;
    global path refHead len LTR dLTR Ts plotCounter
    global Theta
    global Y1 X1
    global v_target
        v_target = 70/3.6;%%%%%%%%%%设置参考车速%%%%%%%%%%%%
    Y1 = 0;
    X1 = 0;
    Theta = zeros(4,1);
    LTR=0;%初始LTR
    dLTR = 0;%初始LTR变化率
    plotCounter = 0;
    path_name = 'new_5pin';
    % path_name = 'largeS';
    % path_name = 'DLC_5pin';
    scale_coeff = 1;
    % 将变量放入基础工作区
    assignin('base', 'path_name', path_name);
    assignin('base', 'scale_coeff', scale_coeff);
    assignin('base','qband_fallback_log',zeros(0,5)); % [t, X, Y, vx, count]
    [path,~,len]=getPath(path_name, scale_coeff);%DLC_trucksim
    figure(1)%************************************************可视化代码
    set(gcf,'Color','white')
    subplot(2,1,1)
    hold on
    plot(path(:,1),path(:,2),'r',LineWidth=1)
    % xlim([1,1300])
    xlabel('X位置 [m]')
    ylabel('Y位置 [m]')
    subplot(2,1,2)
    % xlim([1,1300])
    xlabel('X位置 [m]')
    ylabel('LTR_{est}')
    %获得参考航向角
    diff_x=diff(path(:,1));
    diff_y=diff(path(:,2));
    refHead=zeros(size(path,1),1);
    refHead(1:end-1) = atan2(diff_y,diff_x);
    refHead(end)=refHead(end-1);
    
    Ts = 0.05;
    % Initialize the discrete states.
    str = [];             % Set str to an empty matrix.
    ts  = [Ts 0];       % sample time: [period, offset]
    %End of mdlInitializeSizes
end		 
%% 更新离散状态量
%==============================================================
% Update the discrete states
%==============================================================
function sys = mdlUpdates(t,x,u)
    sys = x;
    %End of mdlUpdate.
end
%% 更新输出量
%==============================================================
% Calculate outputs
%==============================================================
function sys = mdlOutputs(t,x,u)
    global path len refHead U LTR dLTR Ts plotCounter
    global X1 Y1
%% 相关参数定义
%MPC控制器参数
% 观测量：observe=[Y1,X1,psi,vy,r,LTR_est]';6×1  Ny=6
    ux1=abs(u(13));
    q_vy = interp1([0 50 120],[45 38 30],ux1);
    q_r = interp1([0 50 120],[140 180 520],ux1);
    q_ltr = interp1([0 50 120],[4500 8500 12000],ux1);
    Q=diag([6800,60,4700,q_vy,q_r,q_ltr]);
    R=interp1([0 50 120],[3.5e4 1.2e5 3.0e5],ux1);  

    Ny = size(Q,1);
    % IPOPT-NMPC 采用单射击非线性预测；24步在80 km/h下约覆盖27 m。
    Np = 24;
    Nc = 6;
    Nr = ceil(Np*0.85);
    rho=2e6;   %松弛因子系数；保持约束有效，同时避免OSQP数值尺度过大
    incre_MB=3; % 
    incre_CSC=1;%
    

%% 主程序
% TruckSimInput 保持原接口；MPC 只读取 state(1), state(2), state(13:16)。
    state=TruckSimInput(u);
    pos = [state(16),state(15)];
    delta=U;
    % 获得参考量
    [idx,idxend]=findTargetIdx(pos,path,len,abs(state(13)),Ts,Np);
    %MPC控制器
    tic
    [ddelta,observe,epsilon,Y]=MPC_Ctrllor(t,state,delta,path(idx:idxend,1:2),refHead(idx:idxend), ...
        Ts,Q,R,Np,Nc,Nr,rho,incre_MB, incre_CSC,pos);
    iter_time = toc; 
    %更新车辆控制状态
    deltacmd=ddelta+delta;
    deltacmd=sign(deltacmd)*min([abs(deltacmd),15/180*pi]);
    U=deltacmd;
    dLTR = (observe(6)-LTR)/Ts;
    
    LTR=observe(6);
    Y1 = Y(1:Ny:Np*Ny);
    X1 = Y(2:Ny:Np*Ny);

    LTRs = Y(6:6:Np*6);
    plotCounter = plotCounter + 1;
    if mod(plotCounter,5) == 1
        subplot(2,1,1)
        hold on
        plot([pos(1);X1],[pos(2);Y1],'b')
        hold off
        subplot(2,1,2)
        hold on
        plot([pos(1);X1],[LTR;LTRs],'c')
        hold off
    end
    sys=[deltacmd,LTR,epsilon,iter_time]; % steering
% End of mdlOutputs.
end
%% 子函数:MPC控制器
function [ddelta,observe,epsilon,Y]=MPC_Ctrllor(simT,state,delta,refPos,refHead,Ts,Q,R,Np,Nc,Nr,rho,incre_MB, incre_CSC,pos)
    global Theta v_target
    global plotCounter

    persistent uconstrain yconstrain qband nmpcSolver nmpcSolverKey warmNmpc fallbackCount
    if isempty(fallbackCount)
        fallbackCount = 0;
    end
    u=delta;

    rate_delta=35/180*pi*Ts;%限定方向盘转动速度
    uconstrain=[-15*pi/180,15*pi/180, -rate_delta,rate_delta];
    yconstrain=[-0.59,0.59]; % LTR_est soft bound; |LTR_est|=1 表示达到侧翻边界
    Ny = size(Q,1);

    %参考轨迹
    lenRef=size(refPos,1);
    if lenRef < 1
        refPos = repmat(pos,2,1);
        refHead = [state(14);state(14)];
        lenRef = 2;
    elseif lenRef < 2
        refPos = [refPos;refPos];
        refHead = [refHead;refHead];
        lenRef = 2;
    end
    refPosx=interp1(1:lenRef,refPos(:,1),linspace(1,lenRef,Np));
    refPosy=interp1(1:lenRef,refPos(:,2),linspace(1,lenRef,Np));
    if mod(plotCounter,5) == 0
        subplot(2,1,1)
        hold on
        plot(refPosx,refPosy,'g')
        hold off
    end
    yaw0 = state(14);
    % 变成相对的refHead
    refHeads=interp1(1:lenRef,refHead(:,1),linspace(1,lenRef,Np)) - yaw0;
    % 将参考轨迹给旋转一下
    refPathRot = rotatePath([refPosx-pos(1);refPosy-pos(2)]',[0,0],-yaw0)';
    Yr=reshape([refPathRot(2,:);refPathRot(1,:);refHeads;zeros(3,Np)],Ny*Np,1);

    % Q-band ART: tractor-side 2DOF prediction + preview 调权 + yaw/ay residual feedback。
    if isempty(qband) || qband.Np ~= Np || qband.Nc ~= Nc || abs(qband.Ts-Ts) > 1e-12
        qband = qband_default_config(Ts,Np,Nc);
        nmpcSolver = [];
        nmpcSolverKey = '';
        warmNmpc = [];
    end
    plant = tractor2dof_params(state(13));
    qband.model.ayLimit = plant.ayLimit;
    qband = update_qband_runtime(qband,state,delta,refPosx(:),refPosy(:),refHeads(:)+yaw0,plant,Ts);
    [Q,R,yconstrain,qband] = apply_5pin_segment_schedule(Q,R,yconstrain,qband,pos(1));
    [A,B,C,x,observe,ayNow] = Tractor2DOF_Qband_Model(state,delta,Ts,qband,plant);
    Theta = [qband.residual.dy; qband.residual.dr; ayNow; max(abs(qband.wBand))];

    Nu = size(B,2);
    Ny = size(C,1);
    currentKey = sprintf('QbandIPOPT_Np%d_Nc%d_Ts%.6f_nb%d',Np,Nc,Ts,qband.band.nb);
    if isempty(nmpcSolver) || ~strcmp(nmpcSolverKey,currentKey)
        nmpcSolver = build_qband_ipopt_nmpc_solver(Np,Nc,Ts,qband);
        nmpcSolverKey = currentKey;
        warmNmpc = [];
    end

    refAbs = [refPosx(:)'; refPosy(:)'; (refHeads(:)'+yaw0)];
    xNmpc0 = [state(16); state(15); yaw0; state(1); state(2); observe(6); ...
              qband.residual.dy; qband.residual.dr];
    try
        [du,epsilon,Y,warmNmpc] = solve_qband_ipopt_nmpc(nmpcSolver,xNmpc0,u,Q,R, ...
            refAbs,uconstrain,yconstrain,rho,qband,plant,warmNmpc);
    catch ME
        % IPOPT 若在某个采样步失败，保持上一控制量，保证Simulink仿真不中断。
        fallbackCount = fallbackCount + 1;
        try
            fallbackLog = evalin('base','qband_fallback_log');
        catch
            fallbackLog = zeros(0,5);
        end
        fallbackLog(end+1,:) = [simT,pos(1),pos(2),state(13),fallbackCount];
        assignin('base','qband_fallback_log',fallbackLog);
        if mod(fallbackCount,50) == 1
            warning('QbandART:IpoptFallback','IPOPT-NMPC fallback count %d, latest reason: %s',fallbackCount,ME.message);
        end
        du = 0;
        epsilon = 1;
        Y = Yr;
        Y(1:Ny:Np*Ny) = refPosy(:);
        Y(2:Ny:Np*Ny) = refPosx(:);
        Y(3:Ny:Np*Ny) = refAbs(3,:).';
    end

    assignin('base','qband_weights',qband.wBand(:));
    assignin('base','qband_monitor',qband.monitor);
    assignin('base','qband_preview',qband.previewInfo);

    % IPOPT-NMPC 直接输出绝对坐标预测，无需再从车体坐标旋转回全局坐标。
    %获取相对参考量的控制变化量输出
    ddelta=du;
end

%% Q-band ART CasADi controller helpers
function plant = tractor2dof_params(vx)
    vx = max(2.0,abs(vx));

    % Tractor-side equivalent single-track parameters. These are fixed online
    % calibration parameters, not trailer/liquid states.
    plant.m = 1.15e4;
    plant.Iz = 5.3e4;
    plant.lf = 1.385;
    plant.lr = 4.25;
    plant.Cf = 2.2e5;
    plant.Cr = 4.2e5;
    plant.vx = vx;
    plant.g = 9.806;
    plant.ayLimit = 0.32*plant.g;

    plant.a11 = -(plant.Cf+plant.Cr)/(plant.m*vx);
    plant.a12 = (plant.Cr*plant.lr-plant.Cf*plant.lf)/(plant.m*vx)-vx;
    plant.b1 = plant.Cf/plant.m;
    plant.a21 = (plant.Cr*plant.lr-plant.Cf*plant.lf)/(plant.Iz*vx);
    plant.a22 = -(plant.Cf*plant.lf^2+plant.Cr*plant.lr^2)/(plant.Iz*vx);
    plant.b2 = plant.Cf*plant.lf/plant.Iz;
    plant.cAyVy = plant.a11;
    plant.cAyR = plant.a12 + vx;
    plant.cAyDelta = plant.b1;
end

function [A,B,C,x0,observe,ayNow] = Tractor2DOF_Qband_Model(state,delta,Ts,qband,plant)
    vx = max(2.0,abs(state(13)));
    vy = state(1);
    r = state(2);
    dy = qband.residual.dy;
    dr = qband.residual.dr;
    ayNow = plant.cAyVy*vy + plant.cAyR*r + plant.cAyDelta*delta + dy;
    ltrNow = max(-1.2,min(1.2,ayNow/plant.ayLimit));

    % x=[vy; r; psi_rel; Y_rel; X_rel; vx; LTR_est; d_y; d_r].
    % Only vy and r are the 2DOF dynamic states; the remaining states are
    % kinematic integrators, risk proxy, speed hold, and tractor-side residuals.
    Nx = 9;
    A = eye(Nx);
    B = zeros(Nx,1);

    A(1,1) = 1 + Ts*plant.a11;
    A(1,2) = Ts*plant.a12;
    A(1,8) = Ts;
    B(1) = Ts*plant.b1;

    A(2,1) = Ts*plant.a21;
    A(2,2) = 1 + Ts*plant.a22;
    A(2,9) = Ts;
    B(2) = Ts*plant.b2;

    A(3,2) = Ts;
    A(4,1) = Ts;
    A(4,3) = Ts*vx;
    A(5,6) = Ts;
    A(6,6) = 1;

    ltrDecay = qband.model.ltrMemoryDecay;
    A(7,:) = 0;
    A(7,7) = ltrDecay;
    A(7,1) = (1-ltrDecay)*plant.cAyVy/plant.ayLimit;
    A(7,2) = (1-ltrDecay)*plant.cAyR/plant.ayLimit;
    A(7,8) = (1-ltrDecay)/plant.ayLimit;
    B(7) = (1-ltrDecay)*plant.cAyDelta/plant.ayLimit;

    A(8,8) = qband.eso.residualDecay;
    A(9,9) = qband.eso.residualDecay;

    if isfield(qband.residual,'ltrProxy') && isfinite(qband.residual.ltrProxy)
        ltr0 = qband.residual.ltrProxy;
    else
        ltr0 = ltrNow;
    end
    x0 = [vy; r; 0; 0; 0; vx; ltr0; dy; dr];
    C = [0 0 0 1 0 0 0 0 0;  % Y_rel
         0 0 0 0 1 0 0 0 0;  % X_rel
         0 0 1 0 0 0 0 0 0;  % psi_rel
         1 0 0 0 0 0 0 0 0;  % vy
         0 1 0 0 0 0 0 0 0;  % r
         0 0 0 0 0 0 1 0 0]; % LTR_est
    observe = C*x0;
end

function cfg = qband_default_config(Ts,Np,Nc)
    cfg.Ts = Ts;
    cfg.Np = Np;
    cfg.Nc = Nc;

    % 与论文一致的危险频带投影。主晃动约 0.5 Hz，次频约 1.2 Hz，并保留中高频兜底频带。
    cfg.band.edges = [0.35 0.75;
                      0.75 1.45;
                      1.45 2.40;
                      2.40 3.50];
    cfg.band.nb = size(cfg.band.edges,1);
    cfg.band.nFreqPerBand = 5;
    cfg.band.wBase = [1200.0; 600.0; 60.0; 25.0];
    cfg.band.wAdaptMin = [6.0; 2.5; 0.50; 0.30];
    cfg.band.wPreviewTarget = cfg.band.wBase(:);
    cfg.band.wMin = cfg.band.wAdaptMin(:);
    cfg.band.wMax = [1800; 1200; 420; 260];
    cfg.band.alphaUp = 0.75;
    cfg.band.alphaDown = 0.04;
    cfg.band.softEmax = [8.0; 1.0; 0.20; 0.06];
    cfg.band.softSlackWeight = [2.5e4; 1.2e4; 6.0e3; 3.0e3];

    % demo 中有效的 trick：幅值、差分、二阶差分与 Q-band 共同抑制激励。
    cfg.weight.ay = 18.0;
    cfg.weight.day = 36.0;
    cfg.weight.dday = 8.0;
    cfg.ctrl.ayMax = 3.0;

    % Preview 提前调权：由未来参考曲率/横向加速度识别高风险路径段。
    cfg.preview.safeMaxAy = 2.9;
    cfg.preview.safeMeanAbsAy = 0.75;
    cfg.preview.safeMaxDay = 2.8;
    cfg.preview.safeEnergy = [2.5; 0.70; 0.14; 0.04];
    cfg.preview.energyGain = [1.45; 1.2; 0.9; 0.7];
    cfg.preview.accelGain = 0.75;
    cfg.preview.activationSharpness = 2.0;

    % Feedback 频带监测与残余晃动记忆：输入存在且输出被放大的频带会被动态识别为危险频带。
    cfg.monitor.windowSec = 3.0;
    cfg.monitor.epsEnergy = 1e-5;
    cfg.monitor.gBase = [0.06; 0.05; 0.04; 0.03];
    cfg.monitor.gainResidual = [1.0; 0.6; 0.3; 0.2];
    cfg.monitor.memoryDecay = exp(-Ts/2.0);
    cfg.monitor.inputMinEnergy = [0.06; 0.04; 0.03; 0.02];
    cfg.monitor.fbSharpness = 1.4;
    cfg.monitor.fbTargetGain = 0.55;

    cfg.wBand = cfg.band.wAdaptMin(:);
    cfg.monitor = init_qband_monitor(cfg);
    cfg.previewInfo = struct('E',zeros(cfg.band.nb,1),'lambda',zeros(cfg.band.nb,1), ...
                             'aMax',0,'aMean',0,'dAMax',0);
    cfg.model.ayLimit = 0.38*9.806;
    cfg.model.ltrMemoryTau = 0.55;
    cfg.model.ltrMemoryDecay = exp(-Ts/cfg.model.ltrMemoryTau);
    cfg.eso.residualDecay = 0.92;
    cfg.eso.alpha = 0.18;
    cfg.residual.dy = 0;
    cfg.residual.dr = 0;
    cfg.residual.ltrProxy = NaN;
    cfg.residual.yawResidual = 0;
    cfg.residual.ayResidual = 0;
    cfg.lastDelta = 0;
    cfg.lastVy = NaN;
    cfg.lastYawRate = NaN;
end

function mon = init_qband_monitor(cfg)
    nWin = max(8, round(cfg.monitor.windowSec/cfg.Ts));
    mon = cfg.monitor;
    mon.ayBuf = zeros(nWin,1);
    mon.yBuf = zeros(nWin,1);
    mon.ptr = 1;
    mon.filled = 0;
    mon.Eu = zeros(cfg.band.nb,1);
    mon.Ey = zeros(cfg.band.nb,1);
    mon.G = zeros(cfg.band.nb,1);
    mon.coh = zeros(cfg.band.nb,1);
    mon.risk = zeros(cfg.band.nb,1);
    mon.memory = zeros(cfg.band.nb,1);
end

function cfg = update_qband_runtime(cfg,state,delta,refX,refY,refHeadAbs,plant,Ts)
    vx = max(0.5,abs(state(13)));
    yawRate = state(2);
    if isnan(cfg.lastVy)
        ayMon = vx*yawRate;
        yawResidual = 0;
        ayResidual = 0;
    else
        ayMon = (state(1)-cfg.lastVy)/Ts + vx*yawRate;
        vyDotNom = plant.a11*cfg.lastVy + plant.a12*cfg.lastYawRate + plant.b1*cfg.lastDelta + cfg.residual.dy;
        rDotNom = plant.a21*cfg.lastVy + plant.a22*cfg.lastYawRate + plant.b2*cfg.lastDelta + cfg.residual.dr;
        vyDotMeas = (state(1)-cfg.lastVy)/Ts;
        rDotMeas = (yawRate-cfg.lastYawRate)/Ts;
        ayPred = vyDotNom + vx*cfg.lastYawRate;
        yawResidual = rDotMeas - rDotNom;
        ayResidual = ayMon - ayPred;

        cfg.residual.dy = (1-cfg.eso.alpha)*cfg.residual.dy + cfg.eso.alpha*(vyDotMeas - vyDotNom);
        cfg.residual.dr = (1-cfg.eso.alpha)*cfg.residual.dr + cfg.eso.alpha*yawResidual;
        cfg.residual.dy = max(-2.5,min(2.5,cfg.residual.dy));
        cfg.residual.dr = max(-0.8,min(0.8,cfg.residual.dr));
    end
    cfg.lastVy = state(1);
    cfg.lastYawRate = yawRate;
    cfg.lastDelta = delta;
    cfg.residual.yawResidual = yawResidual;
    cfg.residual.ayResidual = ayResidual;
    ltrInstant = max(-1.25,min(1.25,ayMon/plant.ayLimit));
    if isnan(cfg.residual.ltrProxy)
        cfg.residual.ltrProxy = ltrInstant;
    else
        cfg.residual.ltrProxy = cfg.model.ltrMemoryDecay*cfg.residual.ltrProxy + ...
            (1-cfg.model.ltrMemoryDecay)*ltrInstant;
    end

    % Trailer-free feedback monitor: only tractor-side yaw and acceleration residuals are used.
    yMon = yawResidual + 0.35*ayResidual;
    cfg.monitor = update_qband_monitor(cfg.monitor,cfg,ayMon,yMon);

    [wPreview,previewInfo] = preview_qband_weights(cfg,refX,refY,refHeadAbs,vx);
    cfg.previewInfo = previewInfo;

    wRaw = max(cfg.band.wAdaptMin(:),wPreview(:));
    wFb = cfg.band.wAdaptMin(:);
    for ib = 1:cfg.band.nb
        forcedRisk = max(0,cfg.monitor.Eu(ib)/(cfg.monitor.inputMinEnergy(ib)+1e-9)-1) * ...
                     max(0,cfg.monitor.G(ib)/(cfg.monitor.gBase(ib)+1e-9)-1) * ...
                     max(0.15,cfg.monitor.coh(ib));
        residualRisk = sqrt(cfg.monitor.Ey(ib)) * (1-cfg.monitor.coh(ib));
        cfg.monitor.risk(ib) = forcedRisk;
        cfg.monitor.memory(ib) = max(cfg.monitor.memoryDecay*cfg.monitor.memory(ib),forcedRisk);
        fbAct = 1-exp(-cfg.monitor.fbSharpness*cfg.monitor.memory(ib));
        fbTarget = min(cfg.band.wMax(ib),cfg.band.wPreviewTarget(ib)*(1+cfg.monitor.fbTargetGain*cfg.monitor.memory(ib)));
        wFb(ib) = cfg.band.wAdaptMin(ib) + fbAct*(fbTarget-cfg.band.wAdaptMin(ib)) + ...
                  cfg.monitor.gainResidual(ib)*residualRisk;
    end
    wRaw = max(wRaw,wFb(:));
    wRaw = min(max(wRaw,cfg.band.wMin(:)),cfg.band.wMax(:));

    for ib = 1:cfg.band.nb
        if wRaw(ib) > cfg.wBand(ib)
            alpha = cfg.band.alphaUp;
        else
            alpha = cfg.band.alphaDown;
        end
        cfg.wBand(ib) = (1-alpha)*cfg.wBand(ib) + alpha*wRaw(ib);
    end
    cfg.wBand = min(max(cfg.wBand(:),cfg.band.wMin(:)),cfg.band.wMax(:));
end

function mon = update_qband_monitor(mon,cfg,ay,yMon)
    mon.ayBuf(mon.ptr) = ay;
    mon.yBuf(mon.ptr) = yMon;
    mon.ptr = mon.ptr + 1;
    if mon.ptr > numel(mon.ayBuf)
        mon.ptr = 1;
    end
    mon.filled = min(mon.filled+1,numel(mon.ayBuf));
    if mon.filled < max(8,round(0.8/cfg.Ts))
        return;
    end

    if mon.filled < numel(mon.ayBuf)
        a = mon.ayBuf(1:mon.filled);
        y = mon.yBuf(1:mon.filled);
    else
        idx = [mon.ptr:numel(mon.ayBuf),1:mon.ptr-1];
        a = mon.ayBuf(idx);
        y = mon.yBuf(idx);
    end
    a = a(:)-mean(a);
    y = y(:)-mean(y);
    n = numel(a);
    win = 0.5-0.5*cos(2*pi*(0:n-1)'/max(1,n-1));
    a = a.*win;
    y = y.*win;
    tgrid = (0:n-1)'*cfg.Ts;

    for ib = 1:cfg.band.nb
        fList = linspace(cfg.band.edges(ib,1),cfg.band.edges(ib,2),cfg.band.nFreqPerBand);
        Eu = 0; Ey = 0; Cross = 0;
        for jf = 1:numel(fList)
            c = cos(2*pi*fList(jf)*tgrid);
            s = sin(2*pi*fList(jf)*tgrid);
            ac = c.'*a; as = s.'*a;
            yc = c.'*y; ys = s.'*y;
            Eu = Eu + ac^2 + as^2;
            Ey = Ey + yc^2 + ys^2;
            Cross = Cross + ac*yc + as*ys;
        end
        mon.Eu(ib) = Eu/max(1,n^2);
        mon.Ey(ib) = Ey/max(1,n^2);
        mon.G(ib) = sqrt(mon.Ey(ib)/(mon.Eu(ib)+cfg.monitor.epsEnergy));
        mon.coh(ib) = min(1,max(0,Cross^2/(Eu*Ey+cfg.monitor.epsEnergy)));
    end
end

function [wPreview,info] = preview_qband_weights(cfg,refX,refY,refHeadAbs,vx)
    nb = cfg.band.nb;
    kappa = zeros(numel(refHeadAbs),1);
    if numel(refHeadAbs) > 2
        dpsi = local_wrap_to_pi(diff(refHeadAbs(:)));
        ds = hypot(diff(refX(:)),diff(refY(:)));
        ds = max(ds,0.1);
        kappa(1:end-1) = dpsi(:)./ds(:);
        kappa(end) = kappa(end-1);
    end
    aRef = vx^2*kappa(:);
    if numel(aRef) > cfg.Np
        aRef = aRef(1:cfg.Np);
    elseif numel(aRef) < cfg.Np
        aRef(end+1:cfg.Np,1) = aRef(end);
    end
    dARef = diff(aRef);
    aMax = max(abs(aRef));
    aMean = mean(abs(aRef));
    if isempty(dARef)
        dAMax = 0;
    else
        dAMax = max(abs(dARef));
    end
    rAccel = cfg.preview.accelGain * max([max(0,aMax/(cfg.preview.safeMaxAy+1e-9)-1), ...
                                          0.6*max(0,aMean/(cfg.preview.safeMeanAbsAy+1e-9)-1), ...
                                          0.4*max(0,dAMax/(cfg.preview.safeMaxDay+1e-9)-1)]);
    E = qband_energy_numeric(aRef,cfg);
    lambda = zeros(nb,1);
    wPreview = cfg.band.wAdaptMin(:);
    for ib = 1:nb
        rBand = max(0,E(ib)/(cfg.preview.safeEnergy(ib)+1e-9)-1);
        risk = cfg.preview.energyGain(ib)*rBand + rAccel;
        lambda(ib) = 1-exp(-cfg.preview.activationSharpness*max(0,risk));
        wPreview(ib) = cfg.band.wAdaptMin(ib) + lambda(ib)*(cfg.band.wPreviewTarget(ib)-cfg.band.wAdaptMin(ib));
    end
    wPreview = min(max(wPreview(:),cfg.band.wMin(:)),cfg.band.wMax(:));
    info = struct('E',E,'lambda',lambda,'aMax',aMax,'aMean',aMean,'dAMax',dAMax);
end

function E = qband_energy_numeric(aSeq,cfg)
    aSeq = aSeq(:);
    N = numel(aSeq);
    tgrid = (0:N-1)'*cfg.Ts;
    E = zeros(cfg.band.nb,1);
    for ib = 1:cfg.band.nb
        fList = linspace(cfg.band.edges(ib,1),cfg.band.edges(ib,2),cfg.band.nFreqPerBand);
        Ei = 0;
        for jf = 1:numel(fList)
            c = cos(2*pi*fList(jf)*tgrid);
            s = sin(2*pi*fList(jf)*tgrid);
            Ei = Ei + ((c.'*aSeq)^2 + (s.'*aSeq)^2)/max(1,N^2);
        end
        E(ib) = Ei;
    end
end

function qpSolver = build_qband_casadi_qp_solver(nOpt,nCon)
    import casadi.*
    qp = struct();
    qp.h = Sparsity.dense(nOpt,nOpt);
    qp.a = Sparsity.dense(nCon,nOpt);

    opts = struct('printLevel','none');
    qpSolver = conic('QbandART_QP','qpoases',qp,opts);
end

function [du,epsilon,Y] = solve_qband_casadi_qp(qpSolver,a,b,c,x,u,Q,R,Yr, ...
    uconstrain,yconstrain,rho,qband,Tt,Pipi,Nr)
    %% Fast online Q-band ART QP
    % 该求解器沿用原线性预测模型和约束结构，并将 Q-band 投影、ay/dAy/ddAy
    % 抑制项写成二次代价。ay 采用软惩罚而不是硬约束，避免高风险 DLC 段
    % 因瞬时横向加速度过大导致在线 QP 不可行。
    Nx = size(a,1);
    Nu = size(b,2);
    Ny = size(c,1);
    Np = qband.Np;
    Nc = qband.Nc;
    nMove = Nu*Nc;

    Aaug = [a,b;zeros(Nu,Nx),eye(Nu)];
    Baug = [b;eye(Nu)];
    Caug = [c,zeros(Ny,Nu)];
    ksai = [x;u];

    psai = zeros(Ny*Np,Nx+Nu);
    thetaFull = zeros(Ny*Np,Np*Nu);
    for i = 1:Np
        psai((i-1)*Ny+1:i*Ny,:) = Caug*(Aaug^i);
        for j = 1:i
            thetaFull((i-1)*Ny+1:i*Ny,(j-1)*Nu+1:j*Nu) = Caug*(Aaug^(i-j))*Baug;
        end
    end
    theta = thetaFull*Tt;
    E = psai*ksai;

    Qq = build_output_weight_matrix(Q,Np);
    Rr = kron(eye(Nc),R);

    H = [theta'*Qq*theta + Rr, zeros(nMove,1); zeros(1,nMove), min(rho,1e8)];
    g = [theta'*Qq*(E-Yr); 0];

    [ayBase,ayTheta] = build_ay_prediction_qp(Aaug,Baug,ksai,Tt,Np,Nu,Nx,qband);
    [Hq,gq] = qband_qp_terms(ayBase,ayTheta,qband);
    H(1:nMove,1:nMove) = H(1:nMove,1:nMove) + Hq;
    g(1:nMove) = g(1:nMove) + gq;
    H = project_psd_matrix((H+H')/2,1e2);

    AtMove = zeros(Nc);
    for i = 1:Nc
        AtMove(i,1:i) = 1;
    end
    Umin = kron(ones(Nc,1),uconstrain(:,1));
    Umax = kron(ones(Nc,1),uconstrain(:,2));
    dUmin = [kron(ones(Nc,1),uconstrain(:,3));0];
    dUmax = [kron(ones(Nc,1),uconstrain(:,4));1e3];
    Ut = kron(ones(Nc,1),u);

    Ymin = kron(ones(Nr,1),yconstrain(:,1));
    Ymax = kron(ones(Nr,1),yconstrain(:,2));
    AabsU = [kron(AtMove,eye(Nu)),zeros(nMove,1)];
    AltrHi = [Pipi*theta,-ones(Nr,1)];
    AltrLo = [Pipi*theta, ones(Nr,1)];
    At = [AabsU; AltrHi; AltrLo];

    lba = [Umin-Ut; -ones(Nr,1)*1e10; Ymin-Pipi*E];
    uba = [Umax-Ut; Ymax-Pipi*E;  ones(Nr,1)*1e10];

    sol = qpSolver('h',sparse(H),'g',full(g),'a',sparse(At), ...
        'lbx',full(dUmin),'ubx',full(dUmax),'lba',full(lba),'uba',full(uba));
    stats = qpSolver.stats();
    if isfield(stats,'success') && ~stats.success
        if isfield(stats,'return_status')
            error('QbandART:QpStatus','QP status: %s',char(stats.return_status));
        else
            error('QbandART:QpStatus','QP solver returned unsuccessful status.');
        end
    elseif isfield(stats,'return_status')
        status = lower(char(stats.return_status));
        if isempty(strfind(status,'solved')) && isempty(strfind(status,'optimal')) && ...
                isempty(strfind(status,'successful'))
            error('QbandART:QpStatus','QP status: %s',char(stats.return_status));
        end
    end
    z = full(sol.x);
    du = z(1:Nu);
    epsilon = z(end);
    Y = E + theta*z(1:nMove);
    if any(~isfinite(Y)) || numel(Y) ~= numel(Yr)
        Y = Yr(:);
    end
end

function Hpsd = project_psd_matrix(H,floorEig)
    % OSQP 对很小的负特征值很敏感；该投影只修正数值误差导致的非凸判定。
    H = full((H+H')/2);
    if any(~isfinite(H(:)))
        error('QbandART:NonFiniteHessian','Non-finite Hessian entries.');
    end
    [V,D] = eig(H);
    d = real(diag(D));
    if min(d) < floorEig
        Hpsd = V*diag(max(d,floorEig))*V';
        Hpsd = (Hpsd+Hpsd')/2;
    else
        Hpsd = H;
    end
end

function Qq = build_output_weight_matrix(Q,Np)
    % LTR_est 是风险代理和软约束输出，不再作为强制跟踪到 0 的终端目标。
    % 终端增强仅作用在 Y/X/heading 跟踪项上，避免预测末端 LTR 被代价函数人为压低。
    Ny = size(Q,1);
    Qq = zeros(Ny*Np,Ny*Np);
    gamma = 0.97;
    for k = 1:Np
        Qk = Q;
        decay = gamma^(k-1);
        trackIdx = 1:min(3,Ny);
        Qk(trackIdx,trackIdx) = decay*Q(trackIdx,trackIdx);
        if Ny >= 5
            auxIdx = 4:5;
            Qk(auxIdx,auxIdx) = 0.65*decay*Q(auxIdx,auxIdx);
        end
        if Ny >= 6
            Qk(6,6) = 0.02*Q(6,6);
        end
        if k == Np
            Qk(trackIdx,trackIdx) = 8.0*Qk(trackIdx,trackIdx);
        end
        rows = (k-1)*Ny+1:k*Ny;
        Qq(rows,rows) = Qk;
    end
end
function [Q,R,yconstrain,cfg] = apply_5pin_segment_schedule(Q,R,yconstrain,cfg,posX)
    % Segment-aware tuning for the five independent DLCs.
    % 注意：本函数使用"绝对设置"而非"累积乘法"，避免权重随时间持续增大。
    
    persistent initWbandBase initWeightAy initWeightDay initWeightDday
    if isempty(initWbandBase)
        initWbandBase = cfg.band.wBase(:);
        initWeightAy = 18.0;
        initWeightDay = 36.0;
        initWeightDday = 8.0;
    end
    
    xLook = posX + 35;
    if xLook < 325
        trackScale = 2.55;  steerScale = 0.62;  safetyScale = 1.00;  ayScale = 1.00;  yLimit = 0.76;
    elseif xLook < 600
        trackScale = 2.65;  steerScale = 0.62;  safetyScale = 1.15;  ayScale = 1.15;  yLimit = 0.72;
    elseif xLook < 825
        trackScale = 2.12;  steerScale = 0.85;  safetyScale = 1.25;  ayScale = 1.25;  yLimit = 0.70;
    elseif xLook < 1025
        trackScale = 0.95;  steerScale = 1.08;  safetyScale = 2.00;  ayScale = 2.10;  yLimit = 0.62;
    else
        trackScale = 0.35;  steerScale = 2.28;  safetyScale = 2.85;  ayScale = 2.40;  yLimit = 0.48;
    end

    % Q/R 权重：对传入的 Q（本步新计算的）按段调整
    Q(1,1) = Q(1,1)*trackScale;
    Q(2,2) = Q(2,2)*max(0.05,0.75*trackScale);
    Q(3,3) = Q(3,3)*trackScale;
    if size(Q,1) >= 5
        Q(4,4) = Q(4,4)*max(0.7,min(1.3,trackScale));
        Q(5,5) = Q(5,5)*max(0.9,min(1.5,1/sqrt(trackScale)));
    end
    if size(Q,1) >= 6
        Q(6,6) = Q(6,6)*max(2.5,safetyScale);
    end
    R = R*steerScale;
    yconstrain = [-yLimit,yLimit];

    % Q-band 权重：基于初始值 * 段系数（绝对设置，不随时间累积）
    wBandTarget = min(initWbandBase*safetyScale, cfg.band.wMax(:));
    
    % 平滑地向目标值过渡
    alphaBand = 0.15;
    cfg.wBand = (1-alphaBand)*cfg.wBand(:) + alphaBand*wBandTarget;
    cfg.wBand = min(max(cfg.wBand(:),cfg.band.wMin(:)),cfg.band.wMax(:));
    
    % ay/day/dday 权重：同样基于初始值绝对设置
    cfg.weight.ay = initWeightAy*ayScale;
    cfg.weight.day = initWeightDay*sqrt(ayScale);
    cfg.weight.dday = initWeightDday*sqrt(ayScale);
end

function [ayBase,ayTheta] = build_ay_prediction_qp(Aaug,Baug,ksai,Tt,Np,Nu,Nx,cfg)
    nMove = size(Tt,2);
    ayBase = zeros(Np,1);
    ayTheta = zeros(Np,nMove);
    useRiskState = Nx >= 7 && isfield(cfg,'model') && isfield(cfg.model,'ayLimit');

    xPrev = ksai;
    kPrev = zeros(Nx+Nu,nMove);
    for k = 1:Np
        moveMap = Tt((k-1)*Nu+1:k*Nu,:);
        xCurr = Aaug*xPrev;
        kCurr = Aaug*kPrev + Baug*moveMap;

        if useRiskState
            ayBase(k) = cfg.model.ayLimit*xCurr(7);
            ayTheta(k,:) = cfg.model.ayLimit*kCurr(7,:);
        else
            vx = max(0.5,abs(xPrev(13)));
            ayBase(k) = (xCurr(1)-xPrev(1))/cfg.Ts + vx*xPrev(2);
            ayTheta(k,:) = (kCurr(1,:)-kPrev(1,:))/cfg.Ts + vx*kPrev(2,:);
        end

        xPrev = xCurr;
        kPrev = kCurr;
    end
end

function [Hq,gq] = qband_qp_terms(ayBase,ayTheta,cfg)
    nMove = size(ayTheta,2);
    Hq = zeros(nMove,nMove);
    gq = zeros(nMove,1);

    [Hq,gq] = add_sequence_quadratic(Hq,gq,ayTheta,ayBase,cfg.weight.ay);
    if numel(ayBase) >= 2
        D1 = diff(eye(numel(ayBase)),1,1);
        [Hq,gq] = add_sequence_quadratic(Hq,gq,D1*ayTheta,D1*ayBase,cfg.weight.day);
    end
    if numel(ayBase) >= 3
        D2 = diff(eye(numel(ayBase)),2,1);
        [Hq,gq] = add_sequence_quadratic(Hq,gq,D2*ayTheta,D2*ayBase,cfg.weight.dday);
    end

    N = numel(ayBase);
    tgrid = (0:N-1)'*cfg.Ts;
    for ib = 1:cfg.band.nb
        fList = linspace(cfg.band.edges(ib,1),cfg.band.edges(ib,2),cfg.band.nFreqPerBand);
        for jf = 1:numel(fList)
            cVec = cos(2*pi*fList(jf)*tgrid)/max(1,N);
            sVec = sin(2*pi*fList(jf)*tgrid)/max(1,N);
            [Hq,gq] = add_projection_quadratic(Hq,gq,cVec.'*ayTheta,cVec.'*ayBase,cfg.wBand(ib));
            [Hq,gq] = add_projection_quadratic(Hq,gq,sVec.'*ayTheta,sVec.'*ayBase,cfg.wBand(ib));
        end
    end
    Hq = (Hq+Hq')/2;
end

function [H,g] = add_sequence_quadratic(H,g,T,b,w)
    if w <= 0 || isempty(T)
        return;
    end
    H = H + 2*w*(T.'*T);
    g = g + 2*w*(T.'*b);
end

function [H,g] = add_projection_quadratic(H,g,row,offset,w)
    if w <= 0
        return;
    end
    row = full(row(:));
    offset = full(offset);
    H = H + 2*w*(row*row.');
    g = g + 2*w*offset*row;
end

function solver = build_qband_ipopt_nmpc_solver(Np,Nc,Ts,cfg)
    import casadi.*
    opti = Opti();

    U = opti.variable(1,Np);          % absolute steering angle sequence
    Sband = opti.variable(cfg.band.nb,1);
    epsLTR = opti.variable(1,1);

    x0 = opti.parameter(8,1);         % [X; Y; psi; vy; r; LTR_proxy; dy; dr]
    ref = opti.parameter(3,Np);       % [Xref; Yref; psiref]
    plant = opti.parameter(11,1);     % [vx a11 a12 b1 a21 a22 b2 cAyVy cAyR cAyDelta ayLimit]
    uPrev = opti.parameter(1,1);
    qTrack = opti.parameter(6,1);     % [qX qY qPsi qVy qR qLTR]
    rMove = opti.parameter(1,1);
    rho = opti.parameter(1,1);
    wBand = opti.parameter(cfg.band.nb,1);
    eBand = opti.parameter(cfg.band.nb,1);
    wSeq = opti.parameter(3,1);       % [ay dAy ddAy]
    uLimit = opti.parameter(4,1);     % [umin umax dumin dumax]
    yLimit = opti.parameter(2,1);     % LTR lower/upper
    decay = opti.parameter(2,1);      % [LTR memory decay; residual decay]

    xk = x0;
    uLast = uPrev;
    Ycell = cell(Np,1);
    Aycell = cell(Np,1);
    J = 0;

    vx = plant(1);
    a11 = plant(2); a12 = plant(3); b1 = plant(4);
    a21 = plant(5); a22 = plant(6); b2 = plant(7);
    cAyVy = plant(8); cAyR = plant(9); cAyDelta = plant(10);
    ayLimit = plant(11);

    for k = 1:Np
        uk = U(k);
        if k > Nc
            opti.subject_to(uk == U(Nc));
        end
        duk = uk-uLast;
        opti.subject_to(uk >= uLimit(1));
        opti.subject_to(uk <= uLimit(2));
        opti.subject_to(duk >= uLimit(3));
        opti.subject_to(duk <= uLimit(4));

        X = xk(1); Y = xk(2); psi = xk(3);
        vy = xk(4); r = xk(5); ltr = xk(6);
        dy = xk(7); dr = xk(8);

        vyDot = a11*vy + a12*r + b1*uk + dy;
        rDot = a21*vy + a22*r + b2*uk + dr;
        ay = cAyVy*vy + cAyR*r + cAyDelta*uk + dy;
        Xdot = vx*cos(psi) - vy*sin(psi);
        Ydot = vx*sin(psi) + vy*cos(psi);

        xNext = [X + Ts*Xdot;
                 Y + Ts*Ydot;
                 psi + Ts*r;
                 vy + Ts*vyDot;
                 r + Ts*rDot;
                 decay(1)*ltr + (1-decay(1))*ay/ayLimit;
                 decay(2)*dy;
                 decay(2)*dr];

        refX = ref(1,k);
        refY = ref(2,k);
        refPsi = ref(3,k);
        dX = xNext(1)-refX;
        dY = xNext(2)-refY;
        eX = cos(refPsi)*dX + sin(refPsi)*dY;
        eY = -sin(refPsi)*dX + cos(refPsi)*dY;
        ePsiCost = 2*(1-cos(xNext(3)-refPsi));

        gamma = 0.985^(k-1);
        if k == Np
            gamma = 4.0*gamma;
        end
        J = J + gamma*(qTrack(1)*eX^2 + qTrack(2)*eY^2 + qTrack(3)*ePsiCost + ...
                       qTrack(4)*xNext(4)^2 + qTrack(5)*xNext(5)^2 + ...
                       qTrack(6)*0.03*xNext(6)^2);
        J = J + rMove*duk^2 + 0.02*rMove*uk^2;
        J = J + wSeq(1)*ay^2;
        if k >= 2
            dAy = ay - Aycell{k-1};
            J = J + wSeq(2)*dAy^2;
        end
        if k >= 3
            ddAy = ay - 2*Aycell{k-1} + Aycell{k-2};
            J = J + wSeq(3)*ddAy^2;
        end

        opti.subject_to(xNext(6) >= yLimit(1)-epsLTR);
        opti.subject_to(xNext(6) <= yLimit(2)+epsLTR);

        Ycell{k} = [xNext(2); xNext(1); xNext(3); xNext(4); xNext(5); xNext(6)];
        Aycell{k} = ay;
        xk = xNext;
        uLast = uk;
    end

    YPred = vertcat(Ycell{:});
    AyPred = vertcat(Aycell{:});
    opti.subject_to(Sband >= 0);
    opti.subject_to(epsLTR >= 0);
    for ib = 1:cfg.band.nb
        Ei = qband_energy_casadi(AyPred,cfg,ib);
        J = J + wBand(ib)*Ei;
        opti.subject_to(Ei <= eBand(ib) + Sband(ib));
        J = J + cfg.band.softSlackWeight(ib)*Sband(ib)^2;
    end
    J = J + rho*epsLTR^2;
    opti.minimize(J);

    pOpts = struct();
    pOpts.expand = true;
    pOpts.print_time = false;
    sOpts = struct();
    sOpts.max_iter = 45;
    sOpts.print_level = 0;
    sOpts.tol = 1e-4;
    sOpts.acceptable_tol = 2e-3;
    sOpts.acceptable_iter = 5;
    opti.solver('ipopt',pOpts,sOpts);

    solver.opti = opti;
    solver.U = U;
    solver.Sband = Sband;
    solver.epsLTR = epsLTR;
    solver.x0 = x0;
    solver.ref = ref;
    solver.plant = plant;
    solver.uPrev = uPrev;
    solver.qTrack = qTrack;
    solver.rMove = rMove;
    solver.rho = rho;
    solver.wBand = wBand;
    solver.eBand = eBand;
    solver.wSeq = wSeq;
    solver.uLimit = uLimit;
    solver.yLimit = yLimit;
    solver.decay = decay;
    solver.YPred = YPred;
    solver.AyPred = AyPred;
end

function [du,epsilon,Y,warm] = solve_qband_ipopt_nmpc(solver,x0,uPrev,Q,R,refAbs, ...
    uconstrain,yconstrain,rho,qband,plant,warm)
    opti = solver.opti;
    Np = qband.Np;
    plantVec = [plant.vx; plant.a11; plant.a12; plant.b1; plant.a21; plant.a22; ...
                plant.b2; plant.cAyVy; plant.cAyR; plant.cAyDelta; plant.ayLimit];
    qTrack = [max(20,0.35*Q(2,2)); Q(1,1); Q(3,3); Q(4,4); Q(5,5); Q(6,6)];

    opti.set_value(solver.x0,x0(:));
    opti.set_value(solver.ref,refAbs);
    opti.set_value(solver.plant,plantVec);
    opti.set_value(solver.uPrev,uPrev);
    opti.set_value(solver.qTrack,qTrack);
    opti.set_value(solver.rMove,max(10,min(R,2.5e6)));
    opti.set_value(solver.rho,max(1e4,min(rho,5e7)));
    opti.set_value(solver.wBand,qband.wBand(:));
    opti.set_value(solver.eBand,qband.band.softEmax(:));
    opti.set_value(solver.wSeq,[qband.weight.ay; qband.weight.day; qband.weight.dday]);
    opti.set_value(solver.uLimit,[uconstrain(1);uconstrain(2);uconstrain(3);uconstrain(4)]);
    opti.set_value(solver.yLimit,[yconstrain(1);yconstrain(2)]);
    opti.set_value(solver.decay,[qband.model.ltrMemoryDecay;qband.eso.residualDecay]);

    if isempty(warm) || ~isfield(warm,'U') || size(warm.U,2) ~= Np
        opti.set_initial(solver.U,uPrev*ones(1,Np));
        opti.set_initial(solver.Sband,zeros(qband.band.nb,1));
        opti.set_initial(solver.epsLTR,0);
    else
        opti.set_initial(solver.U,[warm.U(:,2:end),warm.U(:,end)]);
        opti.set_initial(solver.Sband,warm.Sband);
        opti.set_initial(solver.epsLTR,warm.epsLTR);
    end

    try
        sol = opti.solve_limited();
        Usol = full(sol.value(solver.U));
        Y = full(sol.value(solver.YPred));
        Ssol = full(sol.value(solver.Sband));
        epsSol = full(sol.value(solver.epsLTR));
    catch ME
        Usol = full(opti.debug.value(solver.U));
        Y = full(opti.debug.value(solver.YPred));
        Ssol = full(opti.debug.value(solver.Sband));
        epsSol = full(opti.debug.value(solver.epsLTR));
        if any(~isfinite(Usol(:)))
            rethrow(ME);
        end
    end
    warm.U = Usol;
    warm.Sband = Ssol;
    warm.epsLTR = epsSol;
    du = Usol(1)-uPrev;
    epsilon = max([warm.epsLTR; warm.Sband(:)]);
    if any(~isfinite(Y)) || numel(Y) ~= 6*Np
        Y = zeros(6*Np,1);
        Y(1:6:end) = refAbs(2,:).';
        Y(2:6:end) = refAbs(1,:).';
        Y(3:6:end) = refAbs(3,:).';
    end
end

function solver = build_qband_casadi_solver(Nx,Nu,Ny,Np,Nc,Ts,cfg)
    import casadi.*
    opti = Opti();
    Uvar = opti.variable(Nu,Np);
    Sband = opti.variable(cfg.band.nb,1);
    epsLTR = opti.variable(1,1);

    Aparam = opti.parameter(Nx,Nx);
    Bparam = opti.parameter(Nx,Nu);
    Cparam = opti.parameter(Ny,Nx);
    x0param = opti.parameter(Nx,1);
    uPrevParam = opti.parameter(Nu,1);
    YrParam = opti.parameter(Ny*Np,1);
    qDiagParam = opti.parameter(Ny,1);
    rParam = opti.parameter(1,1);
    rhoParam = opti.parameter(1,1);
    wBandParam = opti.parameter(cfg.band.nb,1);
    eBandParam = opti.parameter(cfg.band.nb,1);
    uLimitParam = opti.parameter(4,1);
    yLimitParam = opti.parameter(2,1);

    J = 0;
    Ycell = cell(Np,1);
    Aycell = cell(Np,1);
    xk = x0param;
    uLast = uPrevParam;
    for k = 1:Np
        uk = Uvar(:,k);
        if k > Nc
            opti.subject_to(uk == Uvar(:,Nc));
        end
        duk = uk-uLast;
        opti.subject_to(uk >= uLimitParam(1));
        opti.subject_to(uk <= uLimitParam(2));
        opti.subject_to(duk >= uLimitParam(3));
        opti.subject_to(duk <= uLimitParam(4));
        xPrev = xk;
        xk = Aparam*xk + Bparam*uk;
        yk = Cparam*xk;
        Ycell{k} = yk;
        ayk = (xk(1)-xPrev(1))/Ts + xPrev(13)*xPrev(2);
        Aycell{k} = ayk;
        yRefk = YrParam((k-1)*Ny+1:k*Ny);
        gamma = 0.97^(k-1);
        if k == Np
            gamma = gamma*10;
        end
        e = yk-yRefk;
        J = J + gamma*sum1(qDiagParam.*(e.^2)) + rParam*sum1(duk.^2);
        J = J + cfg.weight.ay*ayk^2;
        if k >= 2
            dAy = ayk - Aycell{k-1};
            J = J + cfg.weight.day*dAy^2;
        end
        if k >= 3
            ddAy = ayk - 2*Aycell{k-1} + Aycell{k-2};
            J = J + cfg.weight.dday*ddAy^2;
        end
        opti.subject_to(ayk >= -cfg.ctrl.ayMax);
        opti.subject_to(ayk <=  cfg.ctrl.ayMax);
        opti.subject_to(yk(6) >= yLimitParam(1)-epsLTR);
        opti.subject_to(yk(6) <= yLimitParam(2)+epsLTR);
        uLast = uk;
    end

    YPred = vertcat(Ycell{:});
    AyPred = vertcat(Aycell{:});
    opti.subject_to(Sband >= 0);
    opti.subject_to(epsLTR >= 0);
    for ib = 1:cfg.band.nb
        Ei = qband_energy_casadi(AyPred,cfg,ib);
        J = J + wBandParam(ib)*Ei;
        opti.subject_to(Ei <= eBandParam(ib) + Sband(ib));
        J = J + cfg.band.softSlackWeight(ib)*Sband(ib)^2;
    end
    J = J + rhoParam*epsLTR^2;
    opti.minimize(J);

    pOpts = struct();
    pOpts.expand = true;
    pOpts.print_time = false;
    sOpts = struct();
    sOpts.max_iter = 80;
    sOpts.print_level = 0;
    sOpts.sb = 'yes';
    sOpts.tol = 1e-4;
    sOpts.acceptable_tol = 1e-3;
    opti.solver('ipopt',pOpts,sOpts);

    solver.opti = opti;
    solver.U = Uvar;
    solver.Sband = Sband;
    solver.epsLTR = epsLTR;
    solver.A = Aparam;
    solver.B = Bparam;
    solver.C = Cparam;
    solver.x0 = x0param;
    solver.uPrev = uPrevParam;
    solver.Yr = YrParam;
    solver.Qdiag = qDiagParam;
    solver.R = rParam;
    solver.rho = rhoParam;
    solver.wBand = wBandParam;
    solver.eBand = eBandParam;
    solver.uLimit = uLimitParam;
    solver.yLimit = yLimitParam;
    solver.YPred = YPred;
    solver.AyPred = AyPred;
end

function Ei = qband_energy_casadi(AyPred,cfg,ib)
    import casadi.*
    N = numel(AyPred);
    tgrid = (0:N-1)'*cfg.Ts;
    fList = linspace(cfg.band.edges(ib,1),cfg.band.edges(ib,2),cfg.band.nFreqPerBand);
    Ei = 0;
    for jf = 1:numel(fList)
        c = DM(cos(2*pi*fList(jf)*tgrid));
        s = DM(sin(2*pi*fList(jf)*tgrid));
        ampC = c.'*AyPred;
        ampS = s.'*AyPred;
        Ei = Ei + (ampC^2+ampS^2)/max(1,N^2);
    end
end

function [du,epsilon,Y,warm] = solve_qband_casadi(solver,A,B,C,x,u,Q,R,Yr,uconstrain,yconstrain,rho,qband,warm)
    opti = solver.opti;
    Np = qband.Np;
    opti.set_value(solver.A,A);
    opti.set_value(solver.B,B);
    opti.set_value(solver.C,C);
    opti.set_value(solver.x0,x(:));
    opti.set_value(solver.uPrev,u);
    opti.set_value(solver.Yr,Yr(:));
    opti.set_value(solver.Qdiag,diag(Q));
    opti.set_value(solver.R,max(1,min(R,3e6)));
    opti.set_value(solver.rho,max(1e6,min(rho,1e10)));
    opti.set_value(solver.wBand,qband.wBand(:));
    opti.set_value(solver.eBand,qband.band.softEmax(:));
    opti.set_value(solver.uLimit,[uconstrain(1);uconstrain(2);uconstrain(3);uconstrain(4)]);
    opti.set_value(solver.yLimit,[yconstrain(1);yconstrain(2)]);

    if isempty(warm)
        opti.set_initial(solver.U,u*ones(1,Np));
        opti.set_initial(solver.Sband,zeros(qband.band.nb,1));
        opti.set_initial(solver.epsLTR,0);
    else
        opti.set_initial(solver.U,[warm.U(:,2:end),warm.U(:,end)]);
        opti.set_initial(solver.Sband,warm.Sband);
        opti.set_initial(solver.epsLTR,warm.epsLTR);
    end

    sol = opti.solve();
    Usol = full(sol.value(solver.U));
    Y = full(sol.value(solver.YPred));
    warm.U = Usol;
    warm.Sband = full(sol.value(solver.Sband));
    warm.epsLTR = full(sol.value(solver.epsLTR));
    du = Usol(1)-u;
    epsilon = max([warm.epsLTR; warm.Sband(:)]);
    if any(~isfinite(Y)) || numel(Y) ~= numel(Yr)
        Y = Yr(:);
    end
end

function y = local_wrap_to_pi(x)
    y = mod(x+pi,2*pi)-pi;
end
%% 子函数：输入单位转换*******************************
% 状态量：state==[vy1,df1,F1,dF1,vy2,df2,F2,dF2,...
%                 th1,dth1,th2,dth2,vx1,f1,Y1,X1]';16×1 Nx=16
function state=TruckSimInput(u)
    vy1=u(1)/3.6;
    df1=u(2)/180*pi;
    F1=0;
    dF1=0;
    vy2=0;
    df2=0;
    F2=0;
    dF2=0;

    th1=0;
    dth1=0;
    th2=0;
    dth2=0;
    vx1=abs(u(13))/3.6;
    f1=u(14)/180*pi;
    Y1=u(15);
    X1=u(16);
    state=[vy1,df1,F1,dF1,vy2,df2,F2,dF2,...
          th1,dth1,th2,dth2,vx1,f1,Y1,X1]';%16×1 Nx=16
end

%% 子函数：获取参考轨迹最近的点
function [idx,idxend]=findTargetIdx(pos,path,len,spd,Ts,Np)
    dist = sum([(path(:,1)-pos(1)).^2,(path(:,2)-pos(2)).^2],2);
    [~,idx]=min(dist); %找到距离当前位置最近的一个参考轨迹点的序号和距离
    dist=abs(len-(len(idx)+spd*Ts*Np));
    [~,idxend]=min(dist); %找到距离预测时域终止时车辆预测位置最近的一个参考轨迹点的序号和距离
    idxend = max(idxend,idx);
    idxend = min(idxend,size(path,1));
end

%% 函数：MPC控制器：MPC_Controllor_qpOASES_Ycons*********************************
%****************************使用说明*****************************************************
% ●输入：
%   a,b,c:  为离散形式的模型A,B,C矩阵
%   Q,R:    最优调节的Q，R矩阵，Q为半正定，R为正定
%   x,u:    状态量x(k)和控制量u(k-1)
%   Np,Nc:  Np和Nc分别为预测时域和控制时域个数
%   uconstrain: 控制量及其变化量的限制，形式如下：
%               [u1min u1max du1min du1max;
%               u2min u2max du2min du2max];
%   yconstrain：为观测量，即系统输出的限制，可以设计为硬约束或者软约束，这里使用软约束。
%               假设观测量数量为3个，即Ny=3，则使用样例如下：
%               [y1min y1max;
%                y2min y2max;
%                y3min;y3max];
%   rho:    为松弛因子权重，大于0的数字，数值大表示限制松弛因子，即对输出量的约束更硬。
% ●输出:
%   du：    控制量的变化量，Nu x 1
%****************************************************************************************
% 调教MPC时可用以下参数初始化：
% a=rand(3);b=rand(3,2);c=eye(3);x=rand(3,1);u=[0;0];Q=eye(3);R=0.1*eye(2);rho=5;
% Np=50;Nc=3;Yr=zeros(Np*size(c,1),1);uconstrain=[-1 1 -0.1 0.1; -2 2 -0.2 0.2];
% yconstrain=[-0.1,0.1;-1,1;-1,1]*0.01;

function [du,epsilon,Y]=MPC_Controllor_qpOASES_Ycons(a,b,c,x,u,Q,R,Np,Nc,Yr,uconstrain,yconstrain,rho,Tt,Pipi,Nr)
    %% 模型处理 
    %统计模型状态、控制量和观测量维度
    Nx=size(a,1); %状态量个数
    Nu=size(b,2); %控制量个数
    Ny=size(c,1); %观测量个数
    %构建控制矩阵
    A=[a,b;zeros(Nu,Nx),eye(Nu)]; %(Nx+Nu) x (Nx+Nu)
    B=[b;eye(Nu)];                %(Nx+Nu) x Nu
    C=[c zeros(Ny,Nu)];           %   Ny   x (Nx+Nu)
    %新的控制量为ksai(k)=[x(k),u(k-1)]'
    ksai=[x;u];
    %新的状态空间表达式为：ksai(k+1)=A*ksai(k)+B*du(k)  
    %输出方程为： ita(k)=C*ksai(k)   %Ny x 1
    
    %% 预测输出
    % 获取相关预测矩阵
    psai=zeros(Ny*Np,Nx+Nu); %矩阵psai
    for i=1:Np
        psai(((i-1)*Ny+1):i*Ny,:)=C*A^i;
    end
    theta=zeros(Np*Ny,Np*Nu); %矩阵theta
    for i=1:Np
       for j=1:i
           if j<=Np
           theta(((i-1)*Ny+1):i*Ny,((j-1)*Nu+1):j*Nu)=C*(A^(i-j))*B;
           else
           end
       end
    end
    theta = theta*Tt;
    %输出方程可以写为 Y=psai*ksai(k)+theta*dU  % Ny*Np x 1
    
    %% 控制
    % 变量设置
    E=psai*ksai;
    Qq=build_output_weight_matrix(Q,Np);% 终端惩罚只作用于轨迹跟踪输出
    Rr=kron(eye(Nc),R);
    % 目标函数设计
    % H=theta'*Qq*theta+Rr;
    H=[theta'*Qq*theta+Rr,zeros(Nu*Nc,1);zeros(1,Nu*Nc),rho]; 
    H=(H+H')/2;%保证矩阵对称
    g=[((E-Yr)'*Qq*theta)';0];
    % 约束条件相关矩阵
    At_tmp=zeros(Nc); %下三角方阵
    for i=1:Nc
        At_tmp(i,1:i)=1;
    end
    Nry = 1; % 限制的输出个数
    At=[kron(At_tmp,eye(Nu)),zeros(Nu*Nc,1);
                       Pipi*theta,-ones(Nry*Nr,1);
                       Pipi*theta, ones(Nry*Nr,1)];
    %控制量及其变化量的限制
    Umin=[kron(ones(Nc,1),uconstrain(:,1))];
    Umax=[kron(ones(Nc,1),uconstrain(:,2))];
    dUmin=[kron(ones(Nc,1),uconstrain(:,3));0];
    dUmax=[kron(ones(Nc,1),uconstrain(:,4));1e3];
    %上一时刻的控制量
    Ut=kron(ones(Nc,1),u);
    %输出量约束
    Ymin=kron(ones(Nr,1),yconstrain(:,1));
    Ymax=kron(ones(Nr,1),yconstrain(:,2));
    %开始求解过程
    options = qpOASES_options('default', 'printLevel', 0); 
    % [dU, FVAL, EXITFLAG, iter, lambda] = qpOASES(H, g, At, dUmin, dUmax, Umin-Ut, Umax-Ut, options); %
    [dU, ~, ~, ~, ~] = qpOASES(H, g, At, dUmin, dUmax,...
        [Umin-Ut;ones(Nry*Nr,1)*-1e10;Ymin-Pipi*E], ...
        [Umax-Ut;Ymax-Pipi*E;ones(Nry*Nr,1)*1e10], options);
    du=dU(1:Nu);
    epsilon = dU(end);

    Y = E+theta*dU(1:end-1);
end

%% 函数：SplinePath()****************************
% 通过输入的插值点cv和点数ptNum，根据CV点的距离生成
% 尽可能均匀的插值。
function path = SplinePath(cv,ptNum)
    t0=xy2distance(cv(:,1),cv(:,2));
    t=linspace(t0(1),t0(end),ptNum)';
    path=[pchip(t0,cv(:,1),t),pchip(t0,cv(:,2),t)];
end


%% *****************半挂式液罐车线性模型*******************************
% 模型： 5DOF半挂车+LDP线性双摆=7DOF模型
% 状态量：state==[vy1,df1,F1,dF1,vy2,df2,F2,dF2,...
%                 th1,dth1,th2,dth2,vx1,f1]';14×1 Nx=14
% 观测量：observe=[v1x,f1,dth1,dth2,LTR]';5×1  Ny=5


function [A,B,C]=TrailorTruck_7DOF_LDP(state0, Ts)
    
    % 初始状态定义
    % 定义动力学模型参数
    persistent params
    persistent a b c d e                    %轴距纵向几何参数
    persistent Tw1 Tw2                      %轮距横向几何参数
    persistent h1 h2 h1c h2c                %质心至侧倾轴线距离
                                            %铰接点至侧倾轴线距离
    persistent m1 m2 m1s m2s g;g=9.806;
    persistent I1zz I1xx I1xz I2zz I2xx I2xz  %惯性参量
    persistent k1 k2 k3                       %线性轮胎模型参数
    persistent kr1 kr2 c1 c2 k12              %悬架模型参数
    persistent m0 mp h0 hp lp cd x0 dev x1 x2 %等效摆模型参数
    persistent tau_v
    if isempty(a)
        load('E:\我的文件\大学学习\毕业设计\控制算法设计\液罐车简化模型建立\重新标定车辆模型\best0_ParRecResult_20250309_183708.mat','params')
        ParamsLDP = params(12:end);
        ParamsLDP(17) = 3.5; % dev
        ParamsLDP(18) = 1.7; % cd
        % ParamsLDP = [1.5, 0.13, 0.77, 0.4763, 0.0, 3.5, 0.33];
        %摆模型参数
        g = 9.806;
        exh = params(11); % 调节整体模型的高度
        ParamsLDP(1:2) = ParamsLDP(1:2) + exh;
        [h0,hp,lp,m0,mp,x0,dev,cd] = LDP_getParams(ParamsLDP); 
        x1 = x0 + dev; x2 = x0 - dev;
        %车辆模型参数
        L1=(5+6.27)/2;L2=(10.15+11.5)/2;       %轴距平均
        a=1.385;b=L1-a;c=5.635-a;e=5.5;d=L2-e; %纵向几何
        Tw1=(2.03+1.863*2)/3;Tw2=1.863;        %轮距平均
        hm1s=1.02; %牵引车质心高度
        hm2s=1.00; %50%充液率时半挂车质心高度
        hhitch=1.1;%交接点高度
        
        hroll1=params(9); hroll2=params(10);            %侧倾中心高度
        h1=hm1s-hroll1;h2=hm2s-hroll2;h1c=hhitch-hroll1;h2c=hhitch-hroll2;%高度几何
        m1s=6310;m1=m1s+570+785*2;             %牵引车质量
        % m2s=5925;m2=m2s+665*2;               %半挂车空载质量
        m2s=20387;m2=m2s+665*2;
        I1xx=6879;I1xz=130;I1zz=19665;         %牵引车惯量
        I2xx=9960;I2xz=0;I2zz=179992;        %半挂车惯量
    
        %轮胎模型参数
        k1=-1e4*params(1);k2=-1e4*params(2);k3=-1e4*params(3);             %侧偏刚度
        kr1=params(4)*1e5;kr2=params(5)*1e5;k12=params(6)*1e5; %侧倾刚度
        c1=params(7)*1e3;c2=params(8)*1e3;                    %悬架等效阻尼
        %车辆纵向动力学参数
        tau_v = 0.5;
    end
  
    %% 定义线性动力学模型
    f10 = state0(14);
    v1x = state0(13); v2x = v1x;
    %线性车辆动力学模型
    % M*dX = A0*X + B0*u
  
    m14 = -m1s*h1*c-I1xz;
    m21 = m1*v1x*h1c-m1s*h1*v1x;
    m24 = I1xx+2*m1s*h1^2-m1s*h1*h1c;
    m55 = m2*v2x*h2c-m2s*h2*v2x;
    m58 = I2xx+2*m2s*h2^2-m2s*h2*h2c;
    
    M=[m1*v1x*c, I1zz, 0,   m14,     0,        0,   0,     0;
        m21,    -I1xz, 0,   m24,     0,        0,   0,     0;
       m1*v1x,      0, 0, -m1s*h1, m2*v2x,     0,   0, -m2s*h2;
         0,         0, 0,     0,  m2*v2x*e, -I2zz,  0, I2xz-m2s*h2*e;
         0,         0, 0,     0,    m55,    -I2xz,  0,    m58;
         1,    -c/v1x, 0, -h1c/v1x,  -1,   -e/v2x,  0, h2c/v2x;
         0,         0, 1,     0,     0,        0,   0,     0;
         0,         0, 0,     0,     0,        0,   1,     0   ];
     
     a11 = (c+a)*k1 + (c-b)*k2;
     a12 = a*(c+a)*k1/v1x - b*(c-b)*k2/v1x - m1*v1x*c;
     a22 = (a*k1-b*k2)*h1c/v1x + (m1s*h1-m1*h1c)*v1x;
     a23 = m1s*g*h1 - kr1 -k12;
     a32 = (a*k1-b*k2)/v1x - m1*v1x;
     a36 = -d*k3/v2x - m2*v2x;
     a46 = -d*(e+d)*k3/v2x - m2*v2x*e;
     a56 = (m2s*h2-m2*h2c)*v2x - d*k3*h2c/v2x;
     a57 = m2s*g*h2 - kr2 -k12;
     %  X=   b1     df1    F1  dF1        b2  df2, F2   dF2 
     Acm0=[ a11,    a12,    0,   0,       0,   0,  0,    0;
       (k1+k2)*h1c, a22,  a23, -c1,       0,   0, k12,   0;
         k1+k2,     a32,    0,   0,      k3,  a36, 0,    0;
            0,       0,     0,   0, (e+d)*k3, a46, 0,    0;
            0,       0,    k12,  0,   k3*h2c, a56, a57, -c2;
            0,      -1,     0,   0,       0,   1,  0,    0;
            0,       0,     0,   1,       0,   0,  0,    0;
            0,       0,     0,   0,       0,   0,  0,    1   ];
     
     Bcm0=[-(c+a)*k1, -k1*h1c, -k1, 0,0,0,0,0,0,0,0,0,0,0]';
     % cM1 = 2/(Tw1*(m1 + m2 + m0 + mp));
     % cM2 = 2/(Tw2*(m1 + m2 + m0 + mp));
     % Bcm0=[-(c+a)*k1, -k1*h1c, -k1, 0,0,0,0,0,0,0,0,0,0,0;
     %         1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, cM1,0;
     %        -1, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, cM1,0;
     %         0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, cM2,0;
     %         0, 0, 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, cM2,0;
     %         0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1/tau_v,0]';
     % Bcm0(:,1) = Bcm0(:,1)/c_larger_delta;
     % Bcm0(:,2:5) = Bcm0(:,2:5)/c_smaller_Mz;
       
     %       1    2   3    4    5    6    7    8   9   10   11  12 13 14
    %  X = vy1,  df1,F1, dF1, vy2, df2,  F2, dF2, th1 dth1 th2 dth2 vx f1
    % dX = dvy1;ddf1;dF1;ddF1;dvy2;ddf2;dF2;ddF2;dth1;ddth1;dth2;ddth2;ax;df1 ←M_
    M_ = zeros(14);
    
    M_(9,9) = 1;
    M_(10,:) = [0,0,0,0,v2x/lp,x1/lp,0,hp/lp-2,0,1,0,0,0,0];
    M_(11,11) = 1;
    M_(12,:) = [0,0,0,0,v2x/lp,x2/lp,0,hp/lp-2,0,0,0,1,0,0];
    M_(13,13) = 1;
    M_(14,14) = 1;
    m310 = 1/2*mp*lp;
    mh0p = (m0*h0+mp*hp);
    M_(3,:) = [0,0,0,0,(mp+m0)*v2x,x0*(mp+m0),0,-mh0p,0,m310,0,m310,0,0];
    M_(4,:) = M_(3,:)*e + ...
    [0,0,0,0,-x0*(mp+m0)*v2x,-1/2*(mp+m0)*(x1^2+x2^2),0,x0*mh0p,0,-m310*x1,0,-m310*x2,0,0];
    m510 = -m310*hp;
    M_(5,:) = M_(3,:)*h2c+...
    [0,0,0,0,-mh0p*v2x,-x0*mh0p,0,(m0+h0^2+mp*hp^2),0,m510,0,m510,0,0];

    A_=zeros(14);

    A_(9,10) = 1;
    A_(14,2) = 1;
    A_(10,:) = [0,0,0,0,0,-v2x/lp,-g/lp,0,-g/lp,-cd,0,0,0,0];
    A_(11,12) = 1;
    A_(12,:) = [0,0,0,0,0,-v2x/lp,-g/lp,0,-g/lp,0,0,-cd,0,0];
    % A_(13,13) = -1/tau_v;
    A_(3,:) = [0,0,0,0,0,-(mp+m0)*v2x,0,0,0,0,0,0,0,0];
    A_(4,:) = A_(3,:)*e + x0*(mp + m0)*v2x;
    A_(5,:) = A_(3,:)*h2c+...
    [0,0,0,0,0,mh0p*v2x,mh0p*g,0,-m310*g,0,-m310*g,0,0,0];

     Acm = ([M,zeros(8,6);zeros(6,14)]+M_)\([Acm0,zeros(8,6);zeros(6,14)]+A_);
     Bcm = ([M,zeros(8,6);zeros(6,14)]+M_)\Bcm0;
     
     %% 增广状态量 并转变bi为vyi
     Trans = diag([v1x,1,1,1,v2x,1,1,1,1,1,1,1,1,1]);
     Ac = [Trans*Acm/Trans,zeros(14,2);
           1,0,0,0,0,0,0,0,0,0,0,0,0,v1x,0,0;
           0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;];  %增广状态量后的矩阵

     Bc=[Trans*Bcm;zeros(2,1)];
     %% 离散化线性模型
     [A,B]=c2d_zoh(Ac,Bc,Ts);             %离散化线性模型
     
    
    %% 观测量计算
    %观测量为 Y1 X1 f1 th dth LTR
    C=[zeros(2,14),eye(2);      % Y1 X1 
       zeros(1,13),1,zeros(1,2); % f1
       zeros(1,9),1,zeros(1,6);  % dth1
       zeros(1,11),1,zeros(1,4);  % dth2
      2/(mean([Tw1,Tw2])*(m1+m2)*g)*...
     [0,0,-kr1,-c1,0,0,-kr2,-c2],zeros(1,8)]; % LTR

end
%% 根据变量获取参数
function [h0,hp,lp,m0,mp,x0,dev,cd] = LDP_getParams(Params)
    m = 14462; % TotalWeight
    h0 = Params(1);
    hp = Params(2);
    lp = Params(3);
    m0 = m*Params(4);
    mp = m - m0;
    x0 = Params(5);
    dev = Params(6);
    cd = Params(7);
end

%% 子函数：MB_CSC_get 获得MB和CSC方法所需的压缩映射矩阵
function [Tt,Pipi,Pi] = MB_CSC_get(Np, Nc, Nr, Nu, Ny, incre_MB, incre_CSC)

    Nums = Nums_get(Np, Nc, incre_MB);

    idx=0;
    T = zeros(Np,Nc);
    for i=1:Nc
        T(idx+1:idx+Nums(i),i)=ones(Nums(i),1);
        idx = idx+Nums(i);
    end
    Tt = kron(T,eye(Nu));
    sprintf(['MB压缩映射数组为：' num2str(Nums)])


    Nums = Nums_get(Np, Nr, incre_CSC);

    idx=0; 
    Pi = zeros(Nr,Np);
    for i=1:Nr
        Pi(i,idx+1)=1;
        idx = idx+Nums(i);
    end
    % Pipi = kron(Pi,eye(Ny));
    Pipi = kron(Pi,[0 0 0 0 0 1]);
    sprintf(['CSC压缩映射数组为：' num2str(Nums)])

end


function Nums = Nums_get(Np, Nc, incre_deg)
    % MB_get函数用于生成共Nc个，初值为1，递增，总和为Np的数字序列
    % 输入参数：
    %   - Np: 数字序列的总和
    %   - Nc: 数字序列的长度
    %   - incre_deg: 递增速度的调节因子,越大往后相差越悬殊往前越紧密
    % 输出参数：
    %   - Nums: 生成的数字序列
    Nums = zeros(1, Nc); % 初始化结果数组
    % 将第一个数设置为1
    Nums(1) = 1;
    % 计算递增因子
    total_sum = 0;
    for i = 2:Nc
        total_sum = total_sum + i^incre_deg;
    end
    factor = (Np - Nc) / total_sum; % 减去Nc-1，留出一个单位给每个数字
    
    % 生成递增数字序列
    for i = 2:Nc
        Nums(i) = round(factor * i^incre_deg) + 1; % 加1确保数字大于等于1
    end
    % 调整数字序列，确保总和为Np
    diff = Np - sum(Nums);
    if diff ~= 0
        [~, idx] = max(Nums); % 找到最大值的索引
        Nums(idx) = Nums(idx) + diff; % 调整最大值
    end
end




