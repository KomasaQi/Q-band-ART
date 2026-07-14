%% *****************MPC控制器介绍*********************************
% 模型： 5DOF半挂车+LDP线性化双摆=7DOF模型
% 状态量：state==[vy1,df1,F1,dF1,vy2,df2,F2,dF2,...
%                 th1,dth1,th2,dth2,vx1,f1]';14×1 Nx=14
% 观测量：observe=[Y1,X1,f1,dth1,dth2,LTR]';5×1  Ny=5
% 控制器：对输出进行约束

%% 主函数 根据标志位flag来切换操作具体方程
% MPC_TrajTrackT1_TPSP_LTR_Ctrllor
function [sys,x0,str,ts] = Pin5_MPC_DiffBrk_TrajTrack_Ctrllor_ART(t,x,u,flag)
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
    sizes.NumInputs      = 16;% 状态量：state==[vy1,df1,F1,dF1,vy2,df2,F2,dF2,...
%                                        th1,dth1,th2,dth2,vx1,f1]';14×1
%                                        Nx=14 
%                               还需要额外的Y1 X1来辅助查找参考量
%                               观测量：observe=[v1x,f1,dth1,dth2,LTR]';5×1  Ny=5
    sizes.DirFeedthrough = 1; % Matrix D is non-empty.
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes); 
    x0 =1e-4;   
    global U; % store current ctrl vector: delta_m
    U=0;
    global path refHead len LTR dLTR Ts
    global Theta
    global Y1 X1
    global v_target
        v_target = 70/3.6;%%%%%%%%%%设置参考车速%%%%%%%%%%%%
    Y1 = 0;
    X1 = 0;
    Theta = zeros(4,1);
    LTR=0;%初始LTR
    dLTR = 0;%初始LTR变化率
    path_name = 'new_5pin';
    % path_name = 'largeS';
    % path_name = 'DLC_5pin';
    scale_coeff = 1;
    % 将变量放入基础工作区
    assignin('base', 'path_name', path_name);
    assignin('base', 'scale_coeff', scale_coeff);
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
    ylabel('LTR')
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
    global path len refHead U LTR dLTR Ts
    global X1 Y1
%% 相关参数定义
%MPC控制器参数
% 观测量：observe=[Y1,X1,f1,dth1,dth2,LTR]';5×1  Ny=5
    ux1=u(13);
    q_dth = interp1([0 50 120],[500 500 500]*0.3,ux1);
    Q=diag([4050,50,2550,q_dth,q_dth,150]);
    R=interp1([0 50 120],[1e4 3e6 3e6]*1,ux1);  

    Ny = size(Q,1);
    Np = round(interp1([0 5 50 120],[40 40 40 40],ux1));
    Nc = ceil(Np*0.6);
    Nr = ceil(Np*0.5);
    rho=1e10;   %松弛因子系数
    incre_MB=3; % 
    incre_CSC=1;%
    

%% 主程序
% 状态量：state==[vy1,df1,F1,dF1,  vy2,df2,F2,dF2,...
%                 th,dth,f1,f2,  vx1,vx2,Y1,X1,  Y2,X2]';18×1 Nx=18
    %组合新的状态量：[模型状态量；控制量]’（Nx+Nu）x 1
    state=TruckSimInput(u);
    pos = [state(16),state(15)];
    delta=U;
    % 获得参考量
    [idx,idxend]=findTargetIdx(pos,path,len,state(13),Ts,Np);
    %MPC控制器
    tic
    [ddelta,observe,epsilon,Y]=MPC_Ctrllor(state,delta,path(idx:idxend,1:2),refHead(idx:idxend), ...
        Ts,Q,R,Np,Nc,Nr,rho,incre_MB, incre_CSC,pos);
    iter_time = toc; 
    %更新车辆控制状态
    deltacmd=ddelta+delta;
    deltacmd=sign(deltacmd)*min([abs(deltacmd),15/180*pi]);
    U=deltacmd;
    dLTR = (observe(5)-LTR)/Ts;
    
    LTR=observe(6);
    Y1 = Y(1:Ny:Np*Ny);
    X1 = Y(2:Ny:Np*Ny);

    LTRs = Y(6:6:Np*6);
    subplot(2,1,1)
    hold on
    plot([pos(1);X1],[pos(2);Y1],'b')
    hold off
    subplot(2,1,2)
    hold on
    plot([pos(1);X1],[LTR;LTRs],'c')
    hold off
    sys=[deltacmd,LTR,epsilon,iter_time]; % steering
% End of mdlOutputs.
end
%% 子函数:MPC控制器
function [ddelta,observe,epsilon,Y]=MPC_Ctrllor(state,delta,refPos,refHead,Ts,Q,R,Np,Nc,Nr,rho,incre_MB, incre_CSC,pos)
    global Theta v_target

    % 获得离散化线性模型和观测量
    state0 = [state(1:13);zeros(3,1)];
    [A,B,C] = TrailorTruck_7DOF_LDP(state0,Ts);
    next_state = A*state + B*delta;
    Theta = next_state(9:12);
    observe = C*state;
    persistent Tt Pipi yconstrain uconstrain 
    %MPC控制器相关参数
    x=state0;% 将yaw1变成0，也就是相对值

    u=delta;
    if isempty(yconstrain)
        rate_delta=35/180*pi*Ts;%限定方向盘转动速度
        uconstrain=[-15*pi/180,15*pi/180, -rate_delta,rate_delta];
        % yconstrain=[-1e5,1e5;
        %             -1e5,1e5;
        %             -1e5,1e5;
        %             -25,25;
        %             -25,25;
        %             -0.8,0.8]; % Y1 X1 f1 dth1 dth2 LTR
        yconstrain=[-0.85,0.85]; %  LTR
    end
    Nu = size(B,2);
    Ny = size(C,1);
    [Tt,Pipi] = MB_CSC_get(Np, Nc, Nr, Nu, Ny, incre_MB, incre_CSC);
    %参考轨迹
    lenRef=size(refPos,1);
    refPosx=interp1(1:lenRef,refPos(:,1),linspace(1,lenRef,Np));
    refPosy=interp1(1:lenRef,refPos(:,2),linspace(1,lenRef,Np));
    subplot(2,1,1)
    hold on
    plot(refPosx,refPosy,'g')
    hold off
    yaw0 = state(14);
    % 变成相对的refHead
    refHeads=interp1(1:lenRef,refHead(:,1),linspace(1,lenRef,Np)) - yaw0;
    % 将参考轨迹给旋转一下
    refPathRot = rotatePath([refPosx-pos(1);refPosy-pos(2)]',[0,0],-yaw0)';
    Yr=reshape([refPathRot(2,:);refPathRot(1,:);refHeads;zeros(3,Np)],Ny*Np,1);
    %获得前轮转角变化量
    [du,epsilon,Y]=MPC_Controllor_qpOASES_Ycons(A,B,C,x,u,Q,R,Np,Nc,Yr,uconstrain,yconstrain,rho,Tt,Pipi,Nr);
    % 重新修正Y1 X1
    pathRot = rotatePath([Y(2:Ny:Np*Ny)+pos(1),Y(1:Ny:Np*Ny)+pos(2)],pos,yaw0);
    Y(2:Ny:Np*Ny) = pathRot(:,1); % X1
    Y(1:Ny:Np*Ny) = pathRot(:,2); % Y1
    %获取相对参考量的控制变化量输出
    ddelta=du;  
end
%% 子函数：输入单位转换*******************************
% 状态量：state==[vy1,df1,F1,dF1,vy2,df2,F2,dF2,...
%                 th1,dth1,th2,dth2,vx1,f1,Y1,X1]';16×1 Nx=16
function state=TruckSimInput(u)
    vy1=u(1)/3.6;
    df1=u(2)/180*pi;
    F1=u(3)/180*pi;
    dF1=u(4)/180*pi;
    vy2=u(5)/3.6;
    df2=u(6)/180*pi;
    F2=u(7)/180*pi;
    dF2=u(8)/180*pi;

    coeff_th = 0.1;
    coeff_dth = 1;

    th1=u(9)*coeff_th;
    dth1=u(10)*coeff_dth;
    th2=u(11)*coeff_th;
    dth2=u(12)*coeff_dth;
    vx1=u(13)/3.6;
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
    GAMMA = 0.97; % 预测时域衰减因子
    xNeyeNp=diag(GAMMA.^(0:Np-1));% 加入终端惩罚
    xNeyeNp(end)=xNeyeNp(end)*10;% 加入终端惩罚
    Qq=kron(xNeyeNp,Q);%Qq=kron(eye(Np),Q)
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



