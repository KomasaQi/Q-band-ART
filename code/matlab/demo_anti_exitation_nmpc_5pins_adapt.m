%% demo_anti_exitation_nmpc_5pins_adapt.m
% 基于侧向激励频带抑制的半挂液罐车简化侧向模型 NMPC 示例
%
% 控制对象：
%   - 小车横向运动 + 一阶速度惯性；
%   - 非簧载质量用矩形底盘表示；
%   - 簧上质量/罐体侧倾倒立摆，铰接点固定在非簧载矩形底盘上；
%   - 罐体用随簧上质量侧倾的椭圆表示；
%   - 固定在簧上质量/罐体上的液体等效正单摆；
%   - 真实仿真对象包含 phi 与液体相对摆角 theta 的非线性耦合；
%   - NMPC 内部不显式加入 phi/theta，只对预测 ay 序列做抗激励约束。
%
% 新增功能：
%   1) 基于力矩平衡的 LTR 近似估算与绘制；
%   2) 一键消融实验：比较 ay 幅值约束、ay/差分/频带代价等组合对峰值 LTR 的影响；
%   3) 5组双移线测试：从4s开始，依次执行不同持续时间的双移线，并输出全程轨迹与分段误差；
%   4) 基于预测时域参考轨迹的 preview 自适应频带权重；
%   5) 可选频带能量软约束，用于比较“代价惩罚”和“软安全约束”的差异。
%
% 依赖：CasADi for MATLAB，IPOPT。
%
% 作者：ChatGPT
% 说明：参数为工程演示初值，不代表某一具体车辆的精确标定参数。

clear; clc; close all;

%% -------------------- 0. CasADi 检查 --------------------
try
    import casadi.*
catch ME
    error(['未能导入 CasADi。请先 addpath 到 CasADi MATLAB 路径。\n', ...
           '示例：addpath(''D:/casadi-windows-matlabR2016a-v3.6.5'');\n', ...
           '原始错误：%s'], ME.message);
end

%% -------------------- 1. 参数设置 --------------------
P = default_params();

% ====================== 用户常用开关 ======================
% 单次仿真：false；一键消融实验：true。
% 消融实验会重复构建/求解多个 NMPC，耗时明显增加。
P.exp.runAblation = true;

% 是否在单次仿真结束后播放动画。消融实验默认不播放动画。
P.exp.animateSingleRun = true;

% 可以在这里修改主要参数
% P.traj.t_start = 3.0;
% P.traj.T_lc    = 3.0;
% P.traj.D       = 3.5;
% P.nmpc.N = 50;
% P.sim.dt = 0.05;
% P.ctrl.adaptiveBand = true;

rng(P.sim.rngSeed);

%% -------------------- 2. 仿真/消融入口 --------------------
if P.exp.runAblation
    results = run_ablation_suite(P);
    plot_ablation_results(results, P);
else
    result = run_single_case(P, 'Full anti-excitation NMPC', true);
    plot_results(result, P);
    if P.exp.animateSingleRun
        animate_tanker(result.time, result.X, result.LTR, P);
    end
end

%% ========================================================================
%%                              函数区
%% ========================================================================

function P = default_params()
    g = 9.81;

    % -------------------- 车辆/对象参数 --------------------
    % 参考 29 吨半载半挂液罐车量级，参数为演示用粗略值。
    P.vehicle.g = g;
    P.vehicle.totalMass = 29000;      % kg，演示总质量量级
    P.vehicle.m_u = 4500;             % kg，非簧载/底盘等效质量
    P.vehicle.m_s = 9500;             % kg，空挂/罐体等效簧上质量
    P.vehicle.m_p = 15000;            % kg，半载液体质量

    % 几何参数
    P.vehicle.h_u = 0.75;             % m，侧倾铰点/非簧载高度
    P.vehicle.l_s = 0.85;             % m，侧倾铰点到簧上质心距离

    % 注意：以下几何参数分工不同，避免混用。
    % bodyWidth 只用于非簧载矩形可视化；tankWidth/tankHeight 只用于罐体椭圆可视化；
    % trackWidth 只用于 LTR 分母归一化，不决定动画车身宽度。
    P.vehicle.bodyWidth = 2.55;       % m，车身/非簧载矩形可视化宽度
    P.vehicle.tankWidth = 2.50;       % m，罐体椭圆可视化宽度
    P.vehicle.tankHeight = 1.80;      % m，罐体椭圆可视化高度
    P.vehicle.trackWidth = 1.75;      % m，轮距，仅用于 LTR 力矩归一化

    P.vehicle.h_tankCenter = 1.25;    % m，罐体几何中心到簧上铰接点距离，仅用于可视化
    P.vehicle.h_liqPivot   = 1.55;    % m，液体主晃动中心/摆铰点到簧上铰接点距离
    P.vehicle.h_p = P.vehicle.h_liqPivot; % m，动力学中液体摆铰点高度

    % 液体主晃动频率：默认 0.5 Hz，用单摆长度匹配
    P.vehicle.f_liq = 0.50;           % Hz，主晃动频率
    P.vehicle.omega_liq = 2*pi*P.vehicle.f_liq;
    P.vehicle.l_p = g / P.vehicle.omega_liq^2;  % m，等效摆长
    P.vehicle.zeta_liq = 0.12;        % 液体等效阻尼比，不大但会较快衰减

    % 悬架侧倾目标频率与阻尼
    % 这里通过等效转动惯量和目标侧倾频率反推 k_s, c_s。
    P.vehicle.I_x_body = 10500;       % kg*m^2，空挂/罐体本体绕侧倾轴惯量粗略值
    Ieff_roll = P.vehicle.I_x_body + ...
                P.vehicle.m_s * P.vehicle.l_s^2 + ...
                P.vehicle.m_p * (P.vehicle.h_p^2 + 0.15*P.vehicle.l_p^2);
    P.vehicle.f_roll = 1.10;          % Hz，侧倾模态演示值
    P.vehicle.zeta_roll = 0.25;       % 悬架侧倾阻尼比
    omega_roll = 2*pi*P.vehicle.f_roll;

    % 倒立摆有重力负刚度，因此悬架刚度需包含抵消重力负刚度后的净刚度
    gravDestab = P.vehicle.m_s*g*P.vehicle.l_s + P.vehicle.m_p*g*P.vehicle.h_p;
    kNet = Ieff_roll * omega_roll^2;
    P.vehicle.k_s = kNet + gravDestab;                    % N*m/rad
    P.vehicle.c_s = 2*P.vehicle.zeta_roll*sqrt(Ieff_roll*kNet); % N*m*s/rad

    % 液体铰链阻尼：c_l = 2*zeta*omega*I_p, I_p=m_p*l_p^2
    P.vehicle.c_l = 2 * P.vehicle.zeta_liq * P.vehicle.omega_liq * ...
                    (P.vehicle.m_p * P.vehicle.l_p^2);

    % 速度一阶惯性环节，用于模拟转向到挂车横向加速度建立的滞后
    P.vehicle.tau_v = 0.35;           % s

    % -------------------- LTR 估计参数 --------------------
    P.ltr.enabled = true;
    P.ltr.clipForPlot = 1.25;         % 仅用于图像纵轴/截断显示，原始指标仍保留
    P.ltr.warning = 0.60;
    P.ltr.danger  = 0.85;
    P.ltr.useAbsPeak = true;

    % -------------------- 轨迹参数：5 组双移线 --------------------
    % 从 4 s 开始，执行 5 组不同时长、相同横向位移 3.5 m 的双移线。
    % 每一行：[换过去时间, 保持时间, 换回时间, 组间/组后保持时间]。
    P.traj.mode = 'five_dlc';
    P.traj.D = 3.5;                   % m，每次横向移动距离
    P.traj.t_start = 4.0;             % s，第一组双移线开始时刻；此前 y_ref=0
    P.traj.dlcPhases = [ ...
        5.00, 4.75, 4.50, 5.00;  ... % 双移线1
        4.25, 4.00, 3.75, 4.00;  ... % 双移线2
        3.50, 3.25, 3.00, 3.50;  ... % 双移线3
        2.75, 2.50, 2.25, 2.75;  ... % 双移线4
        2.00, 1.75, 1.50, 5.00];     % 双移线5
    P.traj.segmentInfo = build_five_dlc_segments(P.traj);

    % -------------------- 仿真参数 --------------------
    P.sim.dt = 0.05;                  % s
    P.sim.T_end = P.traj.segmentInfo.totalEndTime; % s，覆盖完整 5 组双移线与最终保持
    P.sim.rngSeed = 1;

    % -------------------- NMPC 参数 --------------------
    P.nmpc.N = 60;                    % 预测步数。dt=0.05,N=60 对应 3 s
    P.nmpc.ipoptMaxIter = 120;
    P.nmpc.ipoptPrint = 0;

    % -------------------- 控制约束与功能开关 --------------------
    P.ctrl.vcmdMin = -3.5;            % m/s
    P.ctrl.vcmdMax =  3.5;            % m/s
    P.ctrl.ayMax   =  3.0;            % m/s^2，横向加速度硬限幅
    P.ctrl.duMax   =  1.2;            % m/s，每个控制步的 v_cmd 变化限幅

    % 消融实验会切换以下开关。
    P.ctrl.useAyHardConstraint = true;
    P.ctrl.useAyCost = true;
    P.ctrl.useDayCost = true;
    P.ctrl.useDdayCost = true;
    P.ctrl.useBand = true;

    % adaptiveBand 为总开关；preview/feedback 分别控制预测型与反馈型权重更新。
    % preview 部分利用预测时域内参考轨迹的 a_y, d a_y 与频带能量提前调权；
    % feedback 部分沿用 roll_rate 频带响应作为兜底修正。
    P.ctrl.adaptiveBand = true;
    P.ctrl.usePreviewAdaptiveBand = true;
    P.ctrl.useFeedbackAdaptiveBand = true;

    % 频带能量软约束开关：若开启，则对 A_y^T Q_i A_y <= E_i,max + s_i 施加软约束。
    P.ctrl.useBandSoftConstraint = false;

    % 代价权重
    P.weight.y      = 120.0;
    P.weight.v      = 6.0;
    P.weight.u      = 0.05;
    P.weight.du     = 2.0;
    P.weight.ay     = 1.8;
    P.weight.day    = 10.0;
    P.weight.dday   = 4.0;
    P.weight.terminalY = 250.0;

    % -------------------- 频带设置 --------------------
    % 根据题设：主晃动 0.5 Hz，次频约 1.2 Hz 且幅值约 15%。
    % 这里将频带分成四段。第一段覆盖主晃动，第二段覆盖次频。
    P.band.edges = [0.35 0.75;   % 主频附近
                    0.75 1.45;   % 次频附近
                    1.45 2.40;   % 中高频耦合
                    2.40 3.50];  % 高频振荡
    P.band.nb = size(P.band.edges,1);
    P.band.nFreqPerBand = 5;

    % 固定 Qband 使用的权重。该组权重在前期可能偏保守，但在 DLC5 高风险段已验证较合适。
    P.band.wBase = [500.0; 200.0; 0.8; 0.5];

    % 自适应 Qband 的低风险初值、目标权重和上限。
    % 低风险使用 wAdaptMin 保证跟踪；preview 识别高风险时向 wPreviewTarget 快速抬升。
    P.band.wAdaptMin = [2.0; 0.8; 0.10; 0.10];
    P.band.wPreviewTarget = P.band.wBase(:);
    P.band.wMin = P.band.wAdaptMin(:);
    P.band.wMax = [650; 330; 180; 120];

    % 非对称权重平滑：风险上升快，风险下降慢，避免“危险后刚结束就立即释放”。
    P.band.alphaUp = 0.75;
    P.band.alphaDown = 0.04;

    % preview 风险参数。safeEnergy 参考 5 组双移线的预测参考加速度频带能量量级设置：
    % DLC4/DLC5 会明显触发，DLC1/DLC2 基本保持低权重，DLC3 轻度触发。
    P.preview.enabled = true;
    P.preview.useBandEnergy = true;
    P.preview.useAccelStats = true;
    P.preview.safeMaxAy = 3.5;          % m/s^2，超过后认为未来参考横向激励偏危险
    P.preview.safeMeanAbsAy = 0.9;      % m/s^2，预测时域内平均激励参考
    P.preview.safeMaxDay = 3.5;         % m/s^2 per step，用离散差分度量参考加速度突变
    P.preview.safeEnergy = [8.0; 1.0; 0.20; 0.06];
    P.preview.energyGain = [1.15; 1.0; 0.75; 0.6];
    P.preview.accelGain = 0.55;
    P.preview.activationSharpness = 2.0;

    % 频带能量软约束阈值。该阈值只约束预测 A_y 序列，不依赖 LTR 真值。
    % 阈值设置为低风险段不活跃、高风险段明显活跃的量级；slack 用于保证优化可行。
    P.band.softEmax = [8.0; 1.0; 0.20; 0.06];
    P.band.softSlackWeight = [2.5e4; 1.2e4; 6.0e3; 3.0e3];

    % 在线频带监测参数
    P.monitor.windowSec = 3.0;
    P.monitor.epsEnergy = 1e-5;
    P.monitor.gBase = [0.06; 0.05; 0.04; 0.03];
    P.monitor.gainForced = [8; 6; 3; 2];
    P.monitor.gainResidual = [1.0; 0.6; 0.3; 0.2];

    % -------------------- 消融实验设置 --------------------
    P.exp.runAblation = false;
    P.exp.animateSingleRun = true;
    P.exp.verboseEachCase = true;
    P.exp.cases = default_ablation_cases();
end

function cases = default_ablation_cases()
    % 消融实验：保留原有单项与组合，并新增 preview adaptive 与 band soft constraint 对比。
    % 字段说明：
    %   adaptiveBand              : 是否启用自适应频带权重总开关
    %   usePreviewAdaptiveBand    : 是否启用基于参考轨迹预测时域的 preview 调权
    %   useFeedbackAdaptiveBand   : 是否启用基于 roll_rate 频带响应的反馈调权
    %   useBandSoftConstraint     : 是否启用频带能量软约束

    c1  = make_case('Baseline: tracking only',          false, false, false, false, false, false, false, false, false);
    c2  = make_case('+ ay hard constraint',             true,  false, false, false, false, false, false, false, false);
    c3  = make_case('+ ay cost only',                   false, true,  false, false, false, false, false, false, false);
    c4  = make_case('+ dAy only',                       false, false, true,  false, false, false, false, false, false);
    c5  = make_case('+ ddAy only',                      false, false, false, true,  false, false, false, false, false);
    c6  = make_case('+ Qband only',                     false, false, false, false, true,  false, false, false, false);
    c7  = make_case('+ ay + dAy',                       true,  true,  true,  false, false, false, false, false, false);
    c8  = make_case('+ ay + dAy + ddAy',                true,  true,  true,  true,  false, false, false, false, false);
    c0  = make_case('+ Qband soft only',                false, false, false, false, true,  false, false, false, true);

    c9  = make_case('+ fixed Qband full',               true,  true,  true,  true,  true,  false, false, false, false);

    % 反馈型 adaptive：保留原有思路，主要用于对比其滞后性。
    c10 = make_case('+ feedback adaptive Qband full',   true,  true,  true,  true,  true,  true,  false, true,  false);

    % 预测型 adaptive：根据预测时域内的参考轨迹加速度/频带能量提前调高 Qband。
    c11 = make_case('+ preview adaptive Qband full',    true,  true,  true,  true,  true,  true,  true,  false, false);

    % 预测 + 反馈：preview 负责提前，feedback 负责真实系统被激发后的兜底修正。
    c12 = make_case('+ preview+feedback adaptive full', true,  true,  true,  true,  true,  true,  true,  true,  false);

    % 频带能量软约束：验证“低风险不强惩罚，高风险触发约束”的效果。
    c13 = make_case('+ fixed Qband + band soft',        true,  true,  true,  true,  true,  false, false, false, true);
    c14 = make_case('+ preview adaptive + band soft',   true,  true,  true,  true,  true,  true,  true,  true,  true);

    cases = [c1, c2, c3, c4, c5, c6, c7, c8, c0, c9, c10, c11, c12, c13, c14];
end

function c = make_case(name, useAyHard, useAy, useDay, useDday, useBand, adaptiveBand, usePreview, useFeedback, useSoft)
    c.name = name;
    c.useAyHardConstraint = useAyHard;
    c.useAyCost = useAy;
    c.useDayCost = useDay;
    c.useDdayCost = useDday;
    c.useBand = useBand;
    c.adaptiveBand = adaptiveBand;
    c.usePreviewAdaptiveBand = usePreview;
    c.useFeedbackAdaptiveBand = useFeedback;
    c.useBandSoftConstraint = useSoft;
end

function results = run_ablation_suite(Pbase)
    cases = Pbase.exp.cases;
    nCase = numel(cases);
    results = cell(nCase,1);

    fprintf('\n========== 开始一键消融实验：%d 个 case ==========' , nCase);
    fprintf('\n说明：每个 case 都重新构建 NMPC 求解器，耗时会明显增加。\n\n');

    for ic = 1:nCase
        P = Pbase;
        P.ctrl.useAyHardConstraint = cases(ic).useAyHardConstraint;
        P.ctrl.useAyCost = cases(ic).useAyCost;
        P.ctrl.useDayCost = cases(ic).useDayCost;
        P.ctrl.useDdayCost = cases(ic).useDdayCost;
        P.ctrl.useBand = cases(ic).useBand;
        P.ctrl.adaptiveBand = cases(ic).adaptiveBand;
        P.ctrl.usePreviewAdaptiveBand = cases(ic).usePreviewAdaptiveBand;
        P.ctrl.useFeedbackAdaptiveBand = cases(ic).useFeedbackAdaptiveBand;
        P.ctrl.useBandSoftConstraint = cases(ic).useBandSoftConstraint;

        fprintf('\n--- Case %d/%d: %s ---\n', ic, nCase, cases(ic).name);
        results{ic} = run_single_case(P, cases(ic).name, Pbase.exp.verboseEachCase);
    end

    fprintf('\n========== 消融实验完成 ==========' );
    print_ablation_table(results);
end

function result = run_single_case(P, caseName, verbose)
    if nargin < 3
        verbose = true;
    end

    fPlant = build_true_plant_dynamics(P);
    solver = build_nmpc_solver(P);

    Tsim = P.sim.T_end;
    dt   = P.sim.dt;
    Nt   = round(Tsim/dt) + 1;
    time = (0:Nt-1)' * dt;

    % 真实对象状态 x = [y; vy; phi; phi_dot; theta; theta_dot]
    %   phi   : 簧上质量/罐体相对竖直向上的侧倾角，phi=0 为直立
    %   theta : 液体相对罐体竖直向下方向的摆角，theta=0 为相对罐体自然下垂
    x = zeros(6,1);
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = deg2rad(0.0);
    x(4) = 0.0;
    x(5) = deg2rad(0.0);
    x(6) = 0.0;

    uPrev = 0.0;

    % 日志
    logX = zeros(Nt,6);
    logU = zeros(Nt,1);
    logAy = zeros(Nt,1);
    logYref = zeros(Nt,1);
    logVref = zeros(Nt,1);
    logWeights = zeros(Nt, P.band.nb);
    logSolveStatus = strings(Nt,1);
    logSolveTime = zeros(Nt,1);
    logLTR = zeros(Nt,1);
    logMroll = zeros(Nt,1);

    % 在线频带权重：固定 Qband 使用 wBase；自适应 Qband 从低风险权重开始。
    if P.ctrl.adaptiveBand && P.ctrl.useBand
        wBand = P.band.wAdaptMin(:);
    else
        wBand = P.band.wBase(:);
    end

    % 频带监测缓冲
    mon = init_band_monitor(P);

    % NMPC warm start
    warm = [];

    if verbose
        fprintf('开始闭环仿真：%s, Nt=%d, dt=%.3f s, N=%d\n', caseName, Nt, dt, P.nmpc.N);
    end

    for k = 1:Nt
        tk = time(k);

        % 参考轨迹预测
        tPred = tk + (0:P.nmpc.N)' * dt;
        [yRefPred, vRefPred] = trajectory_ref(tPred, P.traj);

        % 当前可观测量
        yMeas  = x(1);
        vyMeas = x(2);
        pMeas  = x(4);  % roll_rate 可观测

        % 在求解 NMPC 前，根据预测时域内参考轨迹的未来激励风险提前调整频带权重。
        % 这一步是主动式 preview 调权；反馈监测 mon 只作为兜底修正。
        if P.ctrl.adaptiveBand && P.ctrl.useBand
            wBand = update_band_weights(mon, P, wBand, yRefPred, vRefPred);
        elseif P.ctrl.useBand
            wBand = P.band.wBase(:);
        else
            wBand = P.band.wAdaptMin(:);
        end

        % NMPC 求解
        try
            tic;
            [uCmd, warm, solInfo] = solve_nmpc_once(solver, P, yMeas, vyMeas, uPrev, ...
                                                    yRefPred, vRefPred, wBand, warm);
            logSolveTime(k) = toc;
            logSolveStatus(k) = solInfo.status;
        catch ME
            warning('NMPC 求解失败，使用上一控制量。case=%s, k=%d, t=%.2f, err=%s', ...
                    caseName, k, tk, ME.message);
            uCmd = uPrev;
            logSolveStatus(k) = "failed_use_prev";
            logSolveTime(k) = NaN;
        end

        % 控制限幅，作为最后保护
        uCmd = min(max(uCmd, P.ctrl.vcmdMin), P.ctrl.vcmdMax);

        % 当前实际侧向加速度：由一阶速度环节得到
        ayNow = (uCmd - x(2)) / P.vehicle.tau_v;

        % LTR 估计：基于当前状态与 ay 的侧倾力矩近似
        [ltrNow, mrollNow] = estimate_ltr(x.', ayNow, P);

        % 记录当前时刻
        [yrefNow, vrefNow] = trajectory_ref(tk, P.traj);
        logX(k,:) = x.';
        logU(k) = uCmd;
        logAy(k) = ayNow;
        logYref(k) = yrefNow;
        logVref(k) = vrefNow;
        logWeights(k,:) = wBand(:).';
        logLTR(k) = ltrNow;
        logMroll(k) = mrollNow;

        % 真实对象积分：RK4
        f = @(xx) plant_rhs_numeric(fPlant, xx, uCmd, P);
        xNext = rk4_step(f, x, dt);

        % 更新频带监测器。权重将在下一周期求解前结合 preview 风险一起更新。
        mon = update_band_monitor(mon, P, ayNow, pMeas);

        % 下一步
        x = xNext;
        uPrev = uCmd;

        if verbose && mod(k, max(1,round(1.0/dt))) == 0
            fprintf('t=%.1f / %.1f s, y=%.2f, phi=%.2f deg, theta=%.2f deg, LTR=%.3f, solve=%.3f s\n', ...
                tk, Tsim, x(1), rad2deg(x(3)), rad2deg(x(5)), ltrNow, logSolveTime(k));
        end
    end

    metrics = compute_metrics(time, logX, logU, logAy, logYref, logLTR, logSolveTime, P);

    if verbose
        fprintf('仿真完成：%s。peak |LTR|=%.4f, RMS LTR=%.4f, mean solve=%.4f s\n', ...
            caseName, metrics.peakAbsLTR, metrics.rmsLTR, metrics.meanSolveTime);
    end

    result.name = caseName;
    result.P = P;
    result.time = time;
    result.X = logX;
    result.U = logU;
    result.Ay = logAy;
    result.Yref = logYref;
    result.Vref = logVref;
    result.W = logWeights;
    result.LTR = logLTR;
    result.Mroll = logMroll;
    result.solveStatus = logSolveStatus;
    result.solveTime = logSolveTime;
    result.metrics = metrics;
end

function metrics = compute_metrics(time, X, U, Ay, Yref, LTR, solveTime, P)
    eY = X(:,1) - Yref(:);
    metrics.peakAbsLTR = max(abs(LTR));
    metrics.maxLTR = max(LTR);
    metrics.minLTR = min(LTR);
    metrics.rmsLTR = sqrt(mean(LTR.^2));
    metrics.peakAbsAy = max(abs(Ay));
    metrics.peakAbsPhiDeg = max(abs(rad2deg(X(:,3))));
    metrics.peakAbsThetaRelDeg = max(abs(rad2deg(X(:,5))));
    metrics.rmsTrackError = sqrt(mean(eY.^2));
    metrics.finalTrackError = X(end,1) - Yref(end);
    metrics.totalControlVariation = sum(abs(diff(U)));
    metrics.meanSolveTime = mean(solveTime(~isnan(solveTime)));
    metrics.totalTime = time(end);

    % 每个双移线段单独统计跟踪误差。默认每个 segment 包含：换过去、保持、换回、组间保持。
    if isfield(P, 'traj') && isfield(P.traj, 'segmentInfo')
        seg = P.traj.segmentInfo.segments;
        nSeg = numel(seg);
        metrics.segmentNames = strings(nSeg,1);
        metrics.segmentRmsError = zeros(nSeg,1);
        metrics.segmentPeakAbsError = zeros(nSeg,1);
        metrics.segmentMeanAbsError = zeros(nSeg,1);
        metrics.segmentPeakAbsLTR = zeros(nSeg,1);
        for iseg = 1:nSeg
            idx = (time >= seg(iseg).t0) & (time <= seg(iseg).t1);
            if ~any(idx)
                metrics.segmentNames(iseg) = string(seg(iseg).name);
                metrics.segmentRmsError(iseg) = NaN;
                metrics.segmentPeakAbsError(iseg) = NaN;
                metrics.segmentMeanAbsError(iseg) = NaN;
                metrics.segmentPeakAbsLTR(iseg) = NaN;
            else
                ei = eY(idx);
                metrics.segmentNames(iseg) = string(seg(iseg).name);
                metrics.segmentRmsError(iseg) = sqrt(mean(ei.^2));
                metrics.segmentPeakAbsError(iseg) = max(abs(ei));
                metrics.segmentMeanAbsError(iseg) = mean(abs(ei));
                metrics.segmentPeakAbsLTR(iseg) = max(abs(LTR(idx)));
            end
        end
    else
        metrics.segmentNames = strings(0,1);
        metrics.segmentRmsError = [];
        metrics.segmentPeakAbsError = [];
        metrics.segmentMeanAbsError = [];
        metrics.segmentPeakAbsLTR = [];
    end
end

function print_ablation_table(results)
    n = numel(results);
    names = strings(n,1);
    peakLTR = zeros(n,1);
    rmsLTR = zeros(n,1);
    peakAy = zeros(n,1);
    peakPhi = zeros(n,1);
    rmsTrack = zeros(n,1);
    meanSolve = zeros(n,1);

    for i = 1:n
        names(i) = string(results{i}.name);
        peakLTR(i) = results{i}.metrics.peakAbsLTR;
        rmsLTR(i) = results{i}.metrics.rmsLTR;
        peakAy(i) = results{i}.metrics.peakAbsAy;
        peakPhi(i) = results{i}.metrics.peakAbsPhiDeg;
        rmsTrack(i) = results{i}.metrics.rmsTrackError;
        meanSolve(i) = results{i}.metrics.meanSolveTime;
    end

    T = table(names, peakLTR, rmsLTR, peakAy, peakPhi, rmsTrack, meanSolve, ...
        'VariableNames', {'Case','PeakAbsLTR','RMS_LTR','PeakAbsAy','PeakAbsPhi_deg','RMS_TrackErr_m','MeanSolveTime_s'});
    disp(T);

    if ~isempty(results{1}.metrics.segmentRmsError)
        fprintf('\n========== 每组双移线 RMS 跟踪误差 / 峰值 |LTR| ==========' );
        for iseg = 1:numel(results{1}.metrics.segmentNames)
            segName = results{1}.metrics.segmentNames(iseg);
            segRms = zeros(n,1);
            segLTR = zeros(n,1);
            for i = 1:n
                segRms(i) = results{i}.metrics.segmentRmsError(iseg);
                segLTR(i) = results{i}.metrics.segmentPeakAbsLTR(iseg);
            end
            vRmsName = matlab.lang.makeValidName(['RMS_Error_' char(segName)]);
            vLtrName = matlab.lang.makeValidName(['PeakAbsLTR_' char(segName)]);
            Ts = table(names, segRms, segLTR, ...
                'VariableNames', {'Case', vRmsName, vLtrName});
            fprintf('\n\n%s:\n', char(segName));
            disp(Ts);
        end
    end
end

function fPlant = build_true_plant_dynamics(P)
    import casadi.*

    g   = P.vehicle.g;
    mu  = P.vehicle.m_u;
    ms  = P.vehicle.m_s;
    mp  = P.vehicle.m_p;
    Ix  = P.vehicle.I_x_body;
    hu  = P.vehicle.h_u;
    ls  = P.vehicle.l_s;
    hp  = P.vehicle.h_p;
    lp  = P.vehicle.l_p;
    ks  = P.vehicle.k_s;
    cs  = P.vehicle.c_s;
    cl  = P.vehicle.c_l;

    % 广义坐标 q = [y; phi; theta]
    % theta 为液体相对罐体竖直向下方向的摆角，因此液体绝对摆角 alpha=phi+theta。
    q   = SX.sym('q',3,1);
    dq  = SX.sym('dq',3,1);
    ddq = SX.sym('ddq',3,1);

    y = q(1); phi = q(2); theta = q(3);
    yd = dq(1); phid = dq(2); thetad = dq(3);

    % 位置定义：z 向上；簧上质量位于侧倾铰点上方；液体摆相对罐体向下悬挂
    ys = y + ls*sin(phi);
    zs = hu + ls*cos(phi);

    ypivot = y + hp*sin(phi);
    zpivot = hu + hp*cos(phi);
    alpha = phi + theta;
    yp = ypivot + lp*sin(alpha);
    zp = zpivot - lp*cos(alpha);

    % 速度
    qs = [y; phi; theta];
    dqs = [yd; phid; thetad];
    vs = jacobian([ys; zs], qs) * dqs;
    vp = jacobian([yp; zp], qs) * dqs;

    % 动能与势能
    T = 0.5*mu*yd^2 + 0.5*ms*(vs(1)^2 + vs(2)^2) + 0.5*Ix*phid^2 + ...
        0.5*mp*(vp(1)^2 + vp(2)^2);
    V = ms*g*zs + mp*g*zp + 0.5*ks*phi^2;
    L = T - V;

    dLdq  = jacobian(L, q).';
    dLddq = jacobian(L, dq).';

    ddt_dLddq = jacobian(dLddq, q)*dq + jacobian(dLddq, dq)*ddq;

    % 非保守广义力：悬架侧倾阻尼、液体相对铰链阻尼
    Qnc = SX.zeros(3,1);
    Qnc(2) = -cs*phid;
    Qnc(3) = -cl*thetad;

    R = ddt_dLddq - dLdq - Qnc;

    % 只取 phi/theta 方程，小车 y 加速度由一阶环节给定
    Rpt = R(2:3);

    % Rpt 关于 [ydd; phidd; thetadd] 线性，代入 ydd=ay，解 phidd/thetadd
    ay = SX.sym('ay');
    phidd = SX.sym('phidd');
    thetadd = SX.sym('thetadd');
    ddq_sub = [ay; phidd; thetadd];
    Rsub = substitute(Rpt, ddq, ddq_sub);

    A = jacobian(Rsub, [phidd; thetadd]);
    b = substitute(Rsub, [phidd; thetadd], [0;0]);
    qdd_pt = -solve(A, b);

    fPlant = Function('fPlant', {q, dq, ay}, {qdd_pt}, {'q','dq','ay'}, {'ptdd'});
end

function solver = build_nmpc_solver(P)
    import casadi.*

    N = P.nmpc.N;
    dt = P.sim.dt;
    nb = P.band.nb;

    opti = Opti();

    X = opti.variable(2, N+1);        % [y; vy]
    U = opti.variable(1, N);          % v_y_cmd
    Ay = opti.variable(1, N);         % predicted lateral acceleration, for clearer cost/constraints

    x0      = opti.parameter(2,1);
    uPrevP  = opti.parameter(1,1);
    yRefP   = opti.parameter(N+1,1);
    vRefP   = opti.parameter(N+1,1);
    wBandP  = opti.parameter(nb,1);
    eBandMaxP = opti.parameter(nb,1);  % 频带能量软约束阈值，未启用时仍赋值但不使用

    % 初值约束
    opti.subject_to(X(:,1) == x0);

    % 动力学与约束
    for k = 1:N
        ayk = (U(k) - X(2,k)) / P.vehicle.tau_v;
        xNext = [X(1,k) + dt*X(2,k);
                 X(2,k) + dt*ayk];
        opti.subject_to(X(:,k+1) == xNext);
        opti.subject_to(Ay(k) == ayk);

        opti.subject_to(U(k) >= P.ctrl.vcmdMin);
        opti.subject_to(U(k) <= P.ctrl.vcmdMax);
        if P.ctrl.useAyHardConstraint
            opti.subject_to(Ay(k) >= -P.ctrl.ayMax);
            opti.subject_to(Ay(k) <=  P.ctrl.ayMax);
        end

        if k == 1
            opti.subject_to(U(k)-uPrevP >= -P.ctrl.duMax);
            opti.subject_to(U(k)-uPrevP <=  P.ctrl.duMax);
        else
            opti.subject_to(U(k)-U(k-1) >= -P.ctrl.duMax);
            opti.subject_to(U(k)-U(k-1) <=  P.ctrl.duMax);
        end
    end

    % 代价函数
    J = 0;
    for k = 1:N
        eY = X(1,k) - yRefP(k);
        eV = X(2,k) - vRefP(k);
        J = J + P.weight.y*eY^2 + P.weight.v*eV^2 + P.weight.u*U(k)^2;

        if P.ctrl.useAyCost
            J = J + P.weight.ay*Ay(k)^2;
        end

        if k == 1
            du = U(k) - uPrevP;
        else
            du = U(k) - U(k-1);
        end
        J = J + P.weight.du*du^2;

        if P.ctrl.useDayCost && k >= 2
            dAy = Ay(k) - Ay(k-1);
            J = J + P.weight.day*dAy^2;
        end
        if P.ctrl.useDdayCost && k >= 3
            ddAy = Ay(k) - 2*Ay(k-1) + Ay(k-2);
            J = J + P.weight.dday*ddAy^2;
        end
    end

    % 终端位置误差
    J = J + P.weight.terminalY*(X(1,N+1)-yRefP(N+1))^2;

    % 频带代价/软约束：不增加状态，只对预测 Ay 序列投影到危险频带。
    % 若 useBand=true，则作为二次型代价；若 useBandSoftConstraint=true，则作为软约束。
    bandEnergy = cell(nb,1);
    if P.ctrl.useBand || P.ctrl.useBandSoftConstraint
        tgrid = (0:N-1)' * dt;
        for ib = 1:nb
            Ei = 0;
            fList = linspace(P.band.edges(ib,1), P.band.edges(ib,2), P.band.nFreqPerBand);
            for jf = 1:numel(fList)
                f = fList(jf);
                c = cos(2*pi*f*tgrid);
                s = sin(2*pi*f*tgrid);
                cDM = DM(c(:));
                sDM = DM(s(:));
                ampC = cDM.' * Ay.';
                ampS = sDM.' * Ay.';
                % /N^2 做尺度归一化，避免 N 改变时权重完全失效
                Ei = Ei + (ampC^2 + ampS^2) / (N^2);
            end
            bandEnergy{ib} = Ei;
            if P.ctrl.useBand
                J = J + wBandP(ib) * Ei;
            end
        end
    end

    if P.ctrl.useBandSoftConstraint
        Sband = opti.variable(nb,1);
        opti.subject_to(Sband >= 0);
        for ib = 1:nb
            opti.subject_to(bandEnergy{ib} <= eBandMaxP(ib) + Sband(ib));
            J = J + P.band.softSlackWeight(ib) * Sband(ib)^2;
        end
        solverSband = Sband;
    else
        solverSband = [];
    end

    opti.minimize(J);

    p_opts = struct();
    p_opts.expand = true;
    s_opts = struct();
    s_opts.max_iter = P.nmpc.ipoptMaxIter;
    s_opts.print_level = P.nmpc.ipoptPrint;
    s_opts.sb = 'yes';
    s_opts.tol = 1e-4;
    s_opts.acceptable_tol = 1e-3;
    opti.solver('ipopt', p_opts, s_opts);

    solver.opti = opti;
    solver.X = X;
    solver.U = U;
    solver.Ay = Ay;
    solver.x0 = x0;
    solver.uPrev = uPrevP;
    solver.yRef = yRefP;
    solver.vRef = vRefP;
    solver.wBand = wBandP;
    solver.eBandMax = eBandMaxP;
    solver.Sband = solverSband;
end

function [uCmd, warm, info] = solve_nmpc_once(solver, P, y0, v0, uPrev, yRef, vRef, wBand, warm)
    opti = solver.opti;
    N = P.nmpc.N;

    opti.set_value(solver.x0, [y0; v0]);
    opti.set_value(solver.uPrev, uPrev);
    opti.set_value(solver.yRef, yRef(:));
    opti.set_value(solver.vRef, vRef(:));
    opti.set_value(solver.wBand, wBand(:));
    opti.set_value(solver.eBandMax, P.band.softEmax(:));

    % 初值：若有 warm start，则平移；否则用参考速度初始化
    if isempty(warm)
        Xinit = zeros(2,N+1);
        Xinit(1,:) = yRef(:).';
        Xinit(2,:) = vRef(:).';
        Uinit = vRef(1:N).';
        Ayinit = zeros(1,N);
    else
        Xinit = [warm.X(:,2:end), warm.X(:,end)];
        Uinit = [warm.U(:,2:end), warm.U(:,end)];
        Ayinit = [warm.Ay(:,2:end), warm.Ay(:,end)];
    end

    opti.set_initial(solver.X, Xinit);
    opti.set_initial(solver.U, Uinit);
    opti.set_initial(solver.Ay, Ayinit);

    sol = opti.solve();
    Usol = full(sol.value(solver.U));
    Xsol = full(sol.value(solver.X));
    Aysol = full(sol.value(solver.Ay));

    uCmd = Usol(1);
    warm.X = Xsol;
    warm.U = Usol;
    warm.Ay = Aysol;

    info.status = string(sol.stats.return_status);
end

function dx = plant_rhs_numeric(fPlant, x, uCmd, P)
    vy = x(2);
    phi = x(3);
    phid = x(4);
    theta = x(5);
    thetad = x(6);

    ay = (uCmd - vy) / P.vehicle.tau_v;
    q = [x(1); phi; theta];
    dq = [vy; phid; thetad];
    ptdd = full(fPlant(q, dq, ay));

    dx = zeros(6,1);
    dx(1) = vy;
    dx(2) = ay;
    dx(3) = phid;
    dx(4) = ptdd(1);
    dx(5) = thetad;
    dx(6) = ptdd(2);
end

function xNext = rk4_step(f, x, dt)
    k1 = f(x);
    k2 = f(x + 0.5*dt*k1);
    k3 = f(x + 0.5*dt*k2);
    k4 = f(x + dt*k3);
    xNext = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function info = build_five_dlc_segments(traj)
    % 构造 5 组双移线的时间表。
    phases = traj.dlcPhases;
    nSeg = size(phases,1);
    t = traj.t_start;
    segments = repmat(struct('name','','t0',0,'t1',0,'goStart',0,'goEnd',0, ...
                             'hold1End',0,'backStart',0,'backEnd',0,'hold2End',0), nSeg, 1);
    for i = 1:nSeg
        Tgo   = phases(i,1);
        Thold = phases(i,2);
        Tback = phases(i,3);
        Tgap  = phases(i,4);

        segments(i).name = sprintf('DLC%d', i);
        segments(i).t0 = t;
        segments(i).goStart = t;
        segments(i).goEnd = t + Tgo;
        segments(i).hold1End = segments(i).goEnd + Thold;
        segments(i).backStart = segments(i).hold1End;
        segments(i).backEnd = segments(i).backStart + Tback;
        segments(i).hold2End = segments(i).backEnd + Tgap;
        segments(i).t1 = segments(i).hold2End;
        t = segments(i).t1;
    end
    info.segments = segments;
    info.totalEndTime = segments(end).t1;
end

function [yref, vref] = trajectory_ref(t, traj)
    % 5 组双移线参考轨迹。
    % 每次横向移动均采用五次多项式，保证位置、速度、加速度端点连续为零。
    tInputSize = size(t);
    tv = t(:);
    yref = zeros(size(tv));
    vref = zeros(size(tv));

    if ~isfield(traj, 'mode') || ~strcmpi(traj.mode, 'five_dlc')
        [yref, vref] = quintic_segment_ref(tv, traj.D, traj.T_lc, traj.t_start);
        yref = reshape(yref, tInputSize);
        vref = reshape(vref, tInputSize);
        return;
    end

    D = traj.D;
    seg = traj.segmentInfo.segments;

    % 逐段覆盖。段外默认保持 0；每组结束后回到 0 并保持到下一组。
    for i = 1:numel(seg)
        % 换过去：0 -> D
        idx = (tv >= seg(i).goStart) & (tv <= seg(i).goEnd);
        [yy, vv] = quintic_move(tv(idx), 0, D, seg(i).goStart, seg(i).goEnd);
        yref(idx) = yy; vref(idx) = vv;

        % 第一保持：D
        idx = (tv > seg(i).goEnd) & (tv <= seg(i).hold1End);
        yref(idx) = D; vref(idx) = 0;

        % 换回来：D -> 0
        idx = (tv > seg(i).backStart) & (tv <= seg(i).backEnd);
        [yy, vv] = quintic_move(tv(idx), D, 0, seg(i).backStart, seg(i).backEnd);
        yref(idx) = yy; vref(idx) = vv;

        % 第二保持/组间间隔：0
        idx = (tv > seg(i).backEnd) & (tv <= seg(i).hold2End);
        yref(idx) = 0; vref(idx) = 0;
    end

    % 若预测时域超过最后一个 segment，继续保持 0。
    yref = reshape(yref, tInputSize);
    vref = reshape(vref, tInputSize);
end

function [yref, vref] = quintic_segment_ref(t, D, T, tStart)
    if nargin < 4
        tStart = 0;
    end
    t = t(:);
    tau = t - tStart;
    s = min(max(tau/T, 0), 1);
    yref = D * (10*s.^3 - 15*s.^4 + 6*s.^5);
    dsdt = zeros(size(t));
    active = (tau >= 0) & (tau <= T);
    dsdt(active) = 1/T;
    vref = D * (30*s.^2 - 60*s.^3 + 30*s.^4) .* dsdt;
end

function [y, v] = quintic_move(t, y0, y1, t0, t1)
    t = t(:);
    T = max(t1 - t0, eps);
    s = min(max((t - t0)/T, 0), 1);
    h = 10*s.^3 - 15*s.^4 + 6*s.^5;
    dhds = 30*s.^2 - 60*s.^3 + 30*s.^4;
    y = y0 + (y1-y0)*h;
    v = (y1-y0)*dhds/T;
end

function mon = init_band_monitor(P)
    nWin = max(8, round(P.monitor.windowSec / P.sim.dt));
    mon.nWin = nWin;
    mon.ayBuf = zeros(nWin,1);
    mon.pBuf = zeros(nWin,1);
    mon.ptr = 0;
    mon.filled = 0;
    mon.G = zeros(P.band.nb,1);
    mon.Eu = zeros(P.band.nb,1);
    mon.Ey = zeros(P.band.nb,1);
    mon.coh = zeros(P.band.nb,1);
end

function mon = update_band_monitor(mon, P, ay, p)
    mon.ptr = mod(mon.ptr, mon.nWin) + 1;
    mon.ayBuf(mon.ptr) = ay;
    mon.pBuf(mon.ptr) = p;
    mon.filled = min(mon.filled + 1, mon.nWin);

    if mon.filled < max(8, round(0.8/P.sim.dt))
        return;
    end

    % 按时间顺序取窗口
    if mon.filled < mon.nWin
        ayWin = mon.ayBuf(1:mon.filled);
        pWin  = mon.pBuf(1:mon.filled);
    else
        idx = [mon.ptr+1:mon.nWin, 1:mon.ptr];
        ayWin = mon.ayBuf(idx);
        pWin  = mon.pBuf(idx);
    end

    ayWin = ayWin(:) - mean(ayWin);
    pWin  = pWin(:)  - mean(pWin);
    n = numel(ayWin);
    tt = (0:n-1)' * P.sim.dt;

    % Hann 窗，避免频谱泄漏。这里手写，避免依赖 Signal Processing Toolbox。
    if n > 1
        win = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
    else
        win = 1;
    end
    ayWin = ayWin .* win;
    pWin  = pWin  .* win;

    for ib = 1:P.band.nb
        fList = linspace(P.band.edges(ib,1), P.band.edges(ib,2), P.band.nFreqPerBand);
        Eu = 0; Ey = 0; Cross = 0;
        for jf = 1:numel(fList)
            f = fList(jf);
            c = cos(2*pi*f*tt);
            s = sin(2*pi*f*tt);
            Uc = sum(ayWin .* c);
            Us = sum(ayWin .* s);
            Yc = sum(pWin  .* c);
            Ys = sum(pWin  .* s);
            Eu_j = Uc^2 + Us^2;
            Ey_j = Yc^2 + Ys^2;
            Cross_j = Uc*Yc + Us*Ys;
            Eu = Eu + Eu_j;
            Ey = Ey + Ey_j;
            Cross = Cross + Cross_j;
        end
        mon.Eu(ib) = Eu / max(1,n^2);
        mon.Ey(ib) = Ey / max(1,n^2);
        mon.G(ib) = sqrt(mon.Ey(ib) / (mon.Eu(ib) + P.monitor.epsEnergy));
        mon.coh(ib) = min(1, max(0, Cross^2 / ((Eu*Ey) + P.monitor.epsEnergy)));
    end
end

function wNew = update_band_weights(mon, P, wOld, yRefPred, vRefPred)
    % 频带权重更新：
    %   1) preview 部分：根据预测时域内参考轨迹的 a_y, d a_y 和频带能量提前调权；
    %   2) feedback 部分：根据已观测 roll_rate 的频带响应做兜底修正；
    %   3) 非对称平滑：上升快、下降慢，避免高风险段刚结束就立刻释放。
    nb = P.band.nb;

    % 低风险基础权重
    wRaw = P.band.wAdaptMin(:);

    % -------------------- preview: 未来参考激励风险 --------------------
    if P.ctrl.usePreviewAdaptiveBand && P.preview.enabled
        [wPreview, ~] = preview_band_weights(P, yRefPred, vRefPred);
        wRaw = max(wRaw, wPreview(:));
    end

    % -------------------- feedback: 已观测频带放大风险 --------------------
    if P.ctrl.useFeedbackAdaptiveBand
        wFb = P.band.wAdaptMin(:);
        for ib = 1:nb
            forcedRisk = max(0, mon.G(ib) / (P.monitor.gBase(ib) + 1e-9) - 1) * mon.coh(ib);
            residualRisk = sqrt(mon.Ey(ib)) * (1 - mon.coh(ib));
            wFb(ib) = P.band.wAdaptMin(ib) + ...
                      P.monitor.gainForced(ib)*forcedRisk + ...
                      P.monitor.gainResidual(ib)*residualRisk;
        end
        wRaw = max(wRaw, wFb(:));
    end

    % 限幅
    wRaw = min(max(wRaw, P.band.wMin(:)), P.band.wMax(:));

    % -------------------- 非对称平滑 --------------------
    wOld = wOld(:);
    wNew = zeros(nb,1);
    for ib = 1:nb
        if wRaw(ib) > wOld(ib)
            alpha = P.band.alphaUp;
        else
            alpha = P.band.alphaDown;
        end
        wNew(ib) = (1-alpha)*wOld(ib) + alpha*wRaw(ib);
    end
    wNew = min(max(wNew, P.band.wMin(:)), P.band.wMax(:));
end

function [wPreview, info] = preview_band_weights(P, yRefPred, vRefPred)
    % 根据预测时域内的参考速度序列计算参考横向加速度，并估计未来激励风险。
    % 该函数不依赖真实 LTR，也不依赖当前 roll_rate，主要用于提前在高风险轨迹段前调高 Qband。
    dt = P.sim.dt;
    nb = P.band.nb;
    N = numel(vRefPred) - 1;

    if N <= 2
        wPreview = P.band.wAdaptMin(:);
        info = struct('aMax',0,'aMean',0,'dAMax',0,'E',zeros(nb,1),'lambda',zeros(nb,1));
        return;
    end

    aRef = diff(vRefPred(:)) / dt;       % length N
    dARef = diff(aRef);                  % length N-1

    aMax = max(abs(aRef));
    aMean = mean(abs(aRef));
    if isempty(dARef)
        dAMax = 0;
    else
        dAMax = max(abs(dARef));
    end

    % 参考激励统计风险。只作为全局调权项，不区分频带。
    rAccel = 0;
    if P.preview.useAccelStats
        r1 = max(0, aMax / (P.preview.safeMaxAy + 1e-9) - 1);
        r2 = max(0, aMean / (P.preview.safeMeanAbsAy + 1e-9) - 1);
        r3 = max(0, dAMax / (P.preview.safeMaxDay + 1e-9) - 1);
        rAccel = P.preview.accelGain * max([r1, 0.6*r2, 0.4*r3]);
    end

    % 频带能量风险：哪个频带即将被参考轨迹激发，就优先提高哪个频带权重。
    E = preview_band_energy(P, aRef);
    lambda = zeros(nb,1);
    wPreview = P.band.wAdaptMin(:);

    for ib = 1:nb
        if P.preview.useBandEnergy
            rBand = max(0, E(ib)/(P.preview.safeEnergy(ib)+1e-9) - 1);
        else
            rBand = 0;
        end
        risk = P.preview.energyGain(ib)*rBand + rAccel;

        % 将风险映射到 [0,1]。risk=0 时不激活；risk 越大越接近目标权重。
        lambda(ib) = 1 - exp(-P.preview.activationSharpness * max(0,risk));
        wPreview(ib) = P.band.wAdaptMin(ib) + ...
                       lambda(ib) * (P.band.wPreviewTarget(ib) - P.band.wAdaptMin(ib));
    end

    wPreview = min(max(wPreview(:), P.band.wMin(:)), P.band.wMax(:));

    info.aMax = aMax;
    info.aMean = aMean;
    info.dAMax = dAMax;
    info.E = E;
    info.lambda = lambda;
end

function E = preview_band_energy(P, aSeq)
    % 与 NMPC 内部 Qband 一致的离散正弦/余弦投影能量。
    aSeq = aSeq(:);
    N = numel(aSeq);
    tgrid = (0:N-1)' * P.sim.dt;
    E = zeros(P.band.nb,1);

    for ib = 1:P.band.nb
        fList = linspace(P.band.edges(ib,1), P.band.edges(ib,2), P.band.nFreqPerBand);
        Ei = 0;
        for jf = 1:numel(fList)
            f = fList(jf);
            c = cos(2*pi*f*tgrid);
            s = sin(2*pi*f*tgrid);
            Ei = Ei + ((c.'*aSeq)^2 + (s.'*aSeq)^2) / max(1,N^2);
        end
        E(ib) = Ei;
    end
end

function [LTR, Mroll] = estimate_ltr(X, Ay, P)
    % 基于力矩平衡的 LTR 近似估算。
    %
    % 该估计不是真实轮载测量，而是将当前质心位置、重力与横向惯性力
    % 对地面/轮距中心线产生的侧倾力矩归一化：
    %   LTR_est = M_roll / ((m_total*g*T_track)/2)
    % 其中 M_roll = sum_i m_i*(ay*z_i - g*y_i)。
    % y_i,z_i 均相对底盘中心线/地面参考定义。
    %
    % 物理含义：
    %   |LTR_est| 接近 1 表示仅从静动态力矩平衡角度已接近单侧轮载归零。
    % 注意：该估算未包含轮胎、悬架非线性和真实左右轮载传感器信息。

    phi = X(:,3);
    theta = X(:,5);
    alpha = phi + theta;

    g  = P.vehicle.g;
    mu = P.vehicle.m_u;
    ms = P.vehicle.m_s;
    mp = P.vehicle.m_p;
    hu = P.vehicle.h_u;
    ls = P.vehicle.l_s;
    hp = P.vehicle.h_p;
    lp = P.vehicle.l_p;
    T  = P.vehicle.trackWidth;

    Ay = Ay(:);
    phi = phi(:);
    theta = theta(:); %#ok<NASGU>
    alpha = alpha(:);

    % 各质量相对底盘横向中心线的位置。
    y_u = zeros(size(phi));
    z_u = 0.5*hu*ones(size(phi));

    y_s = ls*sin(phi);
    z_s = hu + ls*cos(phi);

    y_p = hp*sin(phi) + lp*sin(alpha);
    z_p = hu + hp*cos(phi) - lp*cos(alpha);

    Mroll = mu*(Ay.*z_u - g*y_u) + ...
            ms*(Ay.*z_s - g*y_s) + ...
            mp*(Ay.*z_p - g*y_p);

    denom = (mu + ms + mp) * g * T / 2;
    LTR = Mroll ./ max(denom, eps);
end

function plot_results(result, P)
    time = result.time;
    X = result.X;
    U = result.U;
    Ay = result.Ay;
    Yref = result.Yref;
    Vref = result.Vref;
    W = result.W;
    LTR = result.LTR;

    figure('Name','NMPC tracking, anti-excitation and LTR results','Color','w','Position',[60 60 1280 920]);

    subplot(5,2,1);
    plot(time, X(:,1), 'LineWidth', 1.6); hold on;
    plot(time, Yref, '--', 'LineWidth', 1.4);
    grid on; xlabel('t (s)'); ylabel('y (m)'); title('小车横向位置'); legend('actual','ref');

    subplot(5,2,2);
    plot(time, X(:,2), 'LineWidth', 1.6); hold on;
    plot(time, Vref, '--', 'LineWidth', 1.4);
    plot(time, U, ':', 'LineWidth', 1.3);
    grid on; xlabel('t (s)'); ylabel('v_y (m/s)'); title('横向速度与命令'); legend('v actual','v ref','v cmd');

    subplot(5,2,3);
    plot(time, Ay, 'LineWidth', 1.6); hold on;
    if P.ctrl.useAyHardConstraint
        yline(P.ctrl.ayMax, '--'); yline(-P.ctrl.ayMax, '--');
    end
    grid on; xlabel('t (s)'); ylabel('a_y (m/s^2)'); title('实际侧向加速度');

    subplot(5,2,4);
    plot(time, rad2deg(X(:,3)), 'LineWidth', 1.6);
    grid on; xlabel('t (s)'); ylabel('\phi (deg)'); title('簧上质量侧倾角');

    subplot(5,2,5);
    plot(time, rad2deg(X(:,4)), 'LineWidth', 1.6);
    grid on; xlabel('t (s)'); ylabel('p=\dot\phi (deg/s)'); title('侧倾角速度');

    subplot(5,2,6);
    plot(time, rad2deg(X(:,5)), 'LineWidth', 1.6); hold on;
    plot(time, rad2deg(X(:,3)+X(:,5)), '--', 'LineWidth', 1.2);
    grid on; xlabel('t (s)'); ylabel('angle (deg)'); title('液体摆角'); legend('\theta relative','\phi+\theta absolute');

    subplot(5,2,7);
    plot(time, U, 'LineWidth', 1.6);
    grid on; xlabel('t (s)'); ylabel('v_{cmd} (m/s)'); title('NMPC 控制量');

    subplot(5,2,8);
    plot(time, W, 'LineWidth', 1.5);
    grid on; xlabel('t (s)'); ylabel('w_i'); title('在线频带权重');
    legend(compose('Band %d',1:P.band.nb), 'Location','best');

    subplot(5,2,9:10);
    ltrPlot = max(min(LTR, P.ltr.clipForPlot), -P.ltr.clipForPlot);
    plot(time, LTR, 'LineWidth', 1.8); hold on;
    plot(time, ltrPlot, ':', 'LineWidth', 1.0);
    yline(P.ltr.warning, '--', 'LTR warning');
    yline(-P.ltr.warning, '--');
    yline(P.ltr.danger, '-.', 'LTR danger');
    yline(-P.ltr.danger, '-.');
    grid on; xlabel('t (s)'); ylabel('LTR_{est}');
    title(sprintf('估算 LTR：peak |LTR| = %.3f, RMS = %.3f', result.metrics.peakAbsLTR, result.metrics.rmsLTR));
    legend('raw LTR','clipped for visual check','Location','best');
end

function plot_ablation_results(results, P)
    n = numel(results);
    names = strings(n,1);
    peakLTR = zeros(n,1);
    rmsLTR = zeros(n,1);
    rmsTrack = zeros(n,1);
    peakAy = zeros(n,1);

    for i = 1:n
        names(i) = string(results{i}.name);
        peakLTR(i) = results{i}.metrics.peakAbsLTR;
        rmsLTR(i) = results{i}.metrics.rmsLTR;
        rmsTrack(i) = results{i}.metrics.rmsTrackError;
        peakAy(i) = results{i}.metrics.peakAbsAy;
    end

    % 图1：原有统计结果 + LTR 对比。
    figure('Name','Ablation: LTR and tracking summary','Color','w','Position',[70 70 1350 880]);

    subplot(2,2,1);
    bar(peakLTR);
    grid on; ylabel('peak |LTR|'); title('峰值 |LTR| 对比');
    set(gca,'XTick',1:n,'XTickLabel',names,'XTickLabelRotation',30);

    subplot(2,2,2);
    bar(rmsLTR);
    grid on; ylabel('RMS LTR'); title('RMS LTR 对比');
    set(gca,'XTick',1:n,'XTickLabel',names,'XTickLabelRotation',30);

    subplot(2,2,3);
    hold on;
    for i = 1:n
        plot(results{i}.time, results{i}.LTR, 'LineWidth', 1.4);
    end
    yline(P.ltr.warning, '--'); yline(-P.ltr.warning, '--');
    yline(P.ltr.danger, '-.'); yline(-P.ltr.danger, '-.');
    grid on; xlabel('t (s)'); ylabel('LTR_{est}'); title('LTR 时间历程');
    legend(names, 'Location','best');

    subplot(2,2,4);
    yyaxis left;
    bar(peakAy, 0.45);
    ylabel('peak |a_y| (m/s^2)');
    yyaxis right;
    plot(1:n, rmsTrack, 'o-', 'LineWidth', 1.6);
    ylabel('RMS tracking error (m)');
    grid on; title('安全-跟踪折中');
    set(gca,'XTick',1:n,'XTickLabel',names,'XTickLabelRotation',30);

    % 图2：全程轨迹对比。
    figure('Name','Ablation: full trajectory comparison','Color','w','Position',[90 90 1350 760]);
    subplot(2,1,1);
    hold on;
    plot(results{1}.time, results{1}.Yref, 'k--', 'LineWidth', 2.0);
    for i = 1:n
        plot(results{i}.time, results{i}.X(:,1), 'LineWidth', 1.25);
    end
    grid on; xlabel('t (s)'); ylabel('y (m)');
    title('全程横向轨迹对比：5组双移线');
    legend(['reference'; names], 'Location','bestoutside');
    mark_dlc_segments(P, ylim);

    subplot(2,1,2);
    hold on;
    for i = 1:n
        plot(results{i}.time, results{i}.X(:,1)-results{i}.Yref, 'LineWidth', 1.25);
    end
    grid on; xlabel('t (s)'); ylabel('tracking error e_y (m)');
    title('全程跟踪误差对比');
    legend(names, 'Location','bestoutside');
    mark_dlc_segments(P, ylim);

    % 图3：每个双移线分段误差对比。
    if ~isempty(results{1}.metrics.segmentRmsError)
        segNames = results{1}.metrics.segmentNames;
        nSeg = numel(segNames);
        segRms = zeros(nSeg,n);
        segPeak = zeros(nSeg,n);
        segLTR = zeros(nSeg,n);
        for i = 1:n
            segRms(:,i) = results{i}.metrics.segmentRmsError(:);
            segPeak(:,i) = results{i}.metrics.segmentPeakAbsError(:);
            segLTR(:,i) = results{i}.metrics.segmentPeakAbsLTR(:);
        end

        figure('Name','Ablation: segmented tracking error','Color','w','Position',[110 110 1350 820]);
        subplot(3,1,1);
        bar(segRms);
        grid on; ylabel('RMS error (m)'); title('每组双移线 RMS 跟踪误差');
        set(gca,'XTick',1:nSeg,'XTickLabel',segNames);
        legend(names, 'Location','bestoutside');

        subplot(3,1,2);
        bar(segPeak);
        grid on; ylabel('peak |error| (m)'); title('每组双移线峰值跟踪误差');
        set(gca,'XTick',1:nSeg,'XTickLabel',segNames);
        legend(names, 'Location','bestoutside');

        subplot(3,1,3);
        bar(segLTR);
        grid on; ylabel('peak |LTR|'); title('每组双移线峰值 |LTR|');
        set(gca,'XTick',1:nSeg,'XTickLabel',segNames);
        legend(names, 'Location','bestoutside');
    end
end

function mark_dlc_segments(P, ylims)
    if ~isfield(P.traj, 'segmentInfo')
        return;
    end
    seg = P.traj.segmentInfo.segments;
    for i = 1:numel(seg)
        xline(seg(i).t0, ':', sprintf('DLC%d start', i), 'LabelVerticalAlignment','bottom');
        xline(seg(i).t1, ':', sprintf('DLC%d end', i), 'LabelVerticalAlignment','top');
    end
    ylim(ylims);
end

function animate_tanker(time, X, LTR, P)
    fig = figure('Name','Simplified tanker lateral animation','Color','w','Position',[120 120 1050 560]); %#ok<NASGU>
    axis equal; grid on; hold on;
    xlabel('lateral y (m)'); ylabel('height z (m)');
    title('底盘横移 - 簧上侧倾 - 液体单摆晃动');

    yMin = min(X(:,1)) - 2.0;
    yMax = max(X(:,1)) + 2.0;
    zMin = -0.2;
    zMax = P.vehicle.h_u + max(P.vehicle.h_tankCenter + P.vehicle.tankHeight/2, P.vehicle.h_p + P.vehicle.l_p) + 1.0;
    xlim([yMin yMax]); ylim([zMin zMax]);

    % 绘图对象
    hGround = plot([yMin yMax], [0 0], 'k-', 'LineWidth', 1.0); %#ok<NASGU>
    hCartPatch = patch(nan,nan,[0.35 0.35 0.35],'EdgeColor','k','LineWidth',1.2);
    hSpring = plot(nan,nan,'k-','LineWidth',1.2);
    hDamper = plot(nan,nan,'Color',[0.3 0.3 0.3],'LineWidth',2.0);
    hBody   = plot(nan,nan,'-','LineWidth',5,'Color',[0.85 0.25 0.10]);
    hTankPatch = patch(nan,nan,[0.88 0.92 1.00],'EdgeColor',[0.1 0.2 0.5],'LineWidth',2.0,'FaceAlpha',0.35);
    hPend   = plot(nan,nan,'-','LineWidth',2,'Color',[0.1 0.2 0.8]);
    hMass   = plot(nan,nan,'o','MarkerSize',14,'MarkerFaceColor',[0.1 0.4 0.9],'MarkerEdgeColor','k');
    hText   = text(yMin+0.2, zMax-0.3, '', 'FontSize', 11);

    % 动画降采样
    stride = max(1, round(0.05 / mean(diff(time))));

    for k = 1:stride:numel(time)
        y = X(k,1);
        phi = X(k,3);
        theta = X(k,5);
        alpha = phi + theta;

        hu = P.vehicle.h_u;
        ls = P.vehicle.l_s;
        hp = P.vehicle.h_p;
        lp = P.vehicle.l_p;

        cartWidth = P.vehicle.bodyWidth;
        cartHeight = 0.28;
        tankA = P.vehicle.tankWidth/2;
        tankB = P.vehicle.tankHeight/2;

        % 非簧载矩形底盘
        cartY = y;
        cartZ = hu;
        cartX = cartY + [-cartWidth/2, cartWidth/2, cartWidth/2, -cartWidth/2];
        cartZs = cartZ + [-cartHeight/2, -cartHeight/2, cartHeight/2, cartHeight/2];
        set(hCartPatch, 'XData', cartX, 'YData', cartZs);

        % 侧倾铰点在矩形底盘上表面
        pivot0 = [cartY; cartZ + cartHeight/2];
        bodyTop  = pivot0 + [ls*sin(phi); ls*cos(phi)];

        % 罐体椭圆，中心固定在簧上质量上，随 phi 倾斜
        tankCenter = pivot0 + [P.vehicle.h_tankCenter*sin(phi); P.vehicle.h_tankCenter*cos(phi)];
        tankAng = linspace(0, 2*pi, 80);
        R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
        ell = R * [tankA*cos(tankAng); tankB*sin(tankAng)] + tankCenter;

        % 悬架弹簧/阻尼示意：从底盘上表面到车身侧倾轴附近
        springBase = pivot0 + [-0.28; 0];
        springTop  = pivot0 + [-0.28 + 0.20*sin(phi); 0.65*cos(phi)];
        damperBase = pivot0 + [0.28; 0];
        damperTop  = pivot0 + [0.28 + 0.20*sin(phi); 0.65*cos(phi)];
        [sx, sz] = spring_polyline(springBase, springTop, 7, 0.05);

        % 液体摆：铰点固定在簧上质量上，位置与罐体几何中心分开设置，正单摆下垂
        pendPivot = pivot0 + [hp*sin(phi); hp*cos(phi)];
        pendEnd = pendPivot + [lp*sin(alpha); -lp*cos(alpha)];

        set(hSpring, 'XData', sx, 'YData', sz);
        set(hDamper, 'XData', [damperBase(1) damperTop(1)], 'YData', [damperBase(2) damperTop(2)]);
        set(hBody, 'XData', [pivot0(1) bodyTop(1)], 'YData', [pivot0(2) bodyTop(2)]);
        set(hTankPatch, 'XData', ell(1,:), 'YData', ell(2,:));
        set(hPend, 'XData', [pendPivot(1) pendEnd(1)], 'YData', [pendPivot(2) pendEnd(2)]);
        set(hMass, 'XData', pendEnd(1), 'YData', pendEnd(2));
        set(hText, 'String', sprintf('t = %.2f s\ny = %.2f m\nphi = %.2f deg\ntheta = %.2f deg\nLTR = %.3f', ...
            time(k), y, rad2deg(phi), rad2deg(theta), LTR(k)));
        drawnow;
    end
end

function [sx, sz] = spring_polyline(p0, p1, nCoil, amp)
    v = p1 - p0;
    L = norm(v);
    if L < 1e-9
        sx = [p0(1), p1(1)];
        sz = [p0(2), p1(2)];
        return;
    end
    e = v / L;
    n = [-e(2); e(1)];
    nPts = 2*nCoil + 2;
    s = linspace(0,1,nPts);
    offset = zeros(1,nPts);
    offset(2:end-1) = amp * (-1).^(1:nPts-2);
    pts = p0 + v*s + n*offset;
    sx = pts(1,:);
    sz = pts(2,:);
end
