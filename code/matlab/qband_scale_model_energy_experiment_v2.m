%% qband_scale_model_energy_experiment_v2.m
% Qband 频带能量阈值的比例实车模型验证脚本（可读性增强版）
%
% 目的：
%   1) 使用一个与 demo_anti_exitation_nmpc_5pins_adapt 思想一致的简化对象 B；
%   2) 将主晃动频率设置为比例模型更接近的 0.5 Hz；
%   3) 验证：目标频带输入能量小于保守阈值时，动态 LTR 可保持在安全裕度内；
%   4) 同时展示“频带能量”比“总输入能量”更能解释晃动/侧翻风险。
%
% 说明：
%   - 这里不是替代完整多体液罐车模型，而是用于说明 Qband 理论与工程约束含义。
%   - 若你已有从真实 demo 模型线性化得到的 A,B,C,D，可直接替换 build_model_B_scale()。
%   - 所有图中都尽量使用中文标题和阈值线，便于汇报和论文插图筛选。

clear; clc; close all;

%% ======================= 1. 参数区 =======================
P = default_params_scale();

% 构造简化对象 B：ay -> dynamic LTR
[A,B,C,D,info] = build_model_B_scale(P);

% 频率响应与理论阈值
freqHz = linspace(0.05, 2.0, 1600);
Gabs = freq_response_mag(A,B,C,D,2*pi*freqHz);
idxBand = freqHz >= P.bandHz(1) & freqHz <= P.bandHz(2);
[GmaxBand, idxLocal] = max(Gabs(idxBand));
fBand = freqHz(idxBand);
fWorst = fBand(idxLocal);

% 保守阈值：若带内输入 RMS 不超过 a_rms_crit，则动态 LTR 峰值不超过 Mdyn。
% 对单频正弦 ay=A*sin(wt)，输入平均功率 Pavg=A^2/2 = arms^2。
Pcrit = (P.Mdyn / max(GmaxBand, eps))^2;     % mean-square / RMS^2 threshold
Acrit = sqrt(2*Pcrit);                       % sinusoidal amplitude threshold

fprintf('\n========== Qband 频带能量阈值分析 ==========' );
fprintf('\n比例模型主频设置: %.2f Hz', P.f_liq);
fprintf('\n目标频带: [%.2f, %.2f] Hz', P.bandHz(1), P.bandHz(2));
fprintf('\n带内最大 |G(jw)| = %.4f LTR/(m/s^2), 出现在 %.3f Hz', GmaxBand, fWorst);
fprintf('\n允许动态 LTR 裕度 Mdyn = %.3f', P.Mdyn);
fprintf('\n保守输入平均功率阈值 Pcrit = %.4f (m/s^2)^2', Pcrit);
fprintf('\n对应单频正弦幅值阈值 Acrit = %.4f m/s^2\n', Acrit);

%% ======================= 2. 构造可解释测试用例 =======================
T = P.Tsim; dt = P.dt;
t = (0:dt:T).';

% 为了避开初值瞬态影响，输入前 4 s 轻微平滑进入。
ramp = smooth_ramp(t, 3.0, 5.0);

cases = struct([]);

% Case 1：主频附近，低于阈值，理论应安全。
cases(1).name = '安全：0.5Hz 带内输入，低于阈值';
cases(1).short = 'Safe in-band';
cases(1).ay = ramp .* (0.75*Acrit*sin(2*pi*P.f_liq*t));
cases(1).expect = '带内功率低于阈值，动态 LTR 应低于裕度线';

% Case 2：主频附近，高于阈值，容易越界。
cases(2).name = '危险：0.5Hz 带内输入，高于阈值';
cases(2).short = 'Unsafe in-band';
cases(2).ay = ramp .* (1.25*Acrit*sin(2*pi*P.f_liq*t));
cases(2).expect = '带内功率超过阈值，动态 LTR 可能越过裕度线';

% Case 3：多频混合，但带内总功率低于阈值。
A3a = 0.45*Acrit; A3b = 0.35*Acrit; A3c = 0.20*Acrit;
cases(3).name = '安全：混合频率，带内功率受控';
cases(3).short = 'Safe mixed';
cases(3).ay = ramp .* (A3a*sin(2*pi*0.42*t) + A3b*sin(2*pi*0.56*t+0.7) + A3c*sin(2*pi*1.15*t));
cases(3).expect = '虽然包含多个频率，但危险频带能量没有集中超限';

% Case 4：多频混合，带内功率高，演示“输入频率混合”风险。
A4a = 0.75*Acrit; A4b = 0.70*Acrit; A4c = 0.25*Acrit;
cases(4).name = '危险：混合频率，带内能量叠加';
cases(4).short = 'Unsafe mixed';
cases(4).ay = ramp .* (A4a*sin(2*pi*0.43*t) + A4b*sin(2*pi*0.58*t+0.8) + A4c*sin(2*pi*1.20*t+1.5));
cases(4).expect = '单个分量未必特别大，但混合后带内能量叠加';

% Case 5：总功率较大，但远离 0.5Hz；用于说明总能量并不等于危险能量。
% 幅值故意大，但放到 1.3Hz，若对象在该处增益较低，LTR 仍可能较低。
cases(5).name = '对照：总功率较大，但远离主频';
cases(5).short = 'Off-band';
cases(5).ay = ramp .* (1.25*Acrit*sin(2*pi*1.30*t));
cases(5).expect = '总输入功率不低，但不集中在主晃动频带，动态 LTR 风险较低';

%% ======================= 3. 仿真与指标计算 =======================
for i = 1:numel(cases)
    [x,z] = simulate_lti(A,B,C,D,t,cases(i).ay);
    cases(i).x = x;
    cases(i).z = z;
    cases(i).phi = x(:,1);
    cases(i).theta = x(:,3);
    cases(i).peakAbsLTR = max(abs(z(t > P.evalStart)));
    cases(i).rmsLTR = rms_local(z(t > P.evalStart));
    cases(i).totalPower = mean(cases(i).ay(t > P.evalStart).^2);
    cases(i).bandPower = band_power_projection(cases(i).ay, t, P.bandHz, P.nFreqPerBand, P.evalStart);
    cases(i).ratioBand = cases(i).bandPower / Pcrit;
    cases(i).ratioTotal = cases(i).totalPower / Pcrit;
    cases(i).safe = cases(i).peakAbsLTR <= P.Mdyn;
end

print_case_table(cases, Pcrit, P);

%% ======================= 4. 可读性增强图 =======================
plot_frequency_story(freqHz, Gabs, P, fWorst, GmaxBand, Pcrit, Acrit);
plot_case_summary(cases, P);
plot_time_story(t, cases, P);
plot_band_vs_total(cases, Pcrit, P);

%% ========================================================================
%%                              函数区
%% ========================================================================
function P = default_params_scale()
    P.g = 9.81;

    % 这里使用比例模型/实车模型识别意义下的主晃动频率。
    % 用户指出约 0.5 Hz，因此频带围绕 0.5 Hz 设置。
    P.f_liq = 0.50;                  % Hz
    P.zeta_liq = 0.10;

    % 侧倾模态不强行设为主峰；本实验关注 ay -> dynamic LTR 的主危险峰。
    P.f_roll = 1.05;                 % Hz
    P.zeta_roll = 0.28;

    % 耦合项：正值表示液体晃动会向侧倾通道传递能量。
    P.k_couple = 0.38;

    % 从 ay 到动态 LTR 的静态/动态标定系数，主要用于形成可解释量级。
    P.gain_phi_to_ltr   = 0.78;
    P.gain_theta_to_ltr = 0.30;
    P.gain_ay_to_ltr    = 0.035;

    % 目标频带。若频率辨识后峰值略偏，可改为 [0.35 0.70] 或 [0.30 0.80]。
    P.bandHz = [0.35, 0.75];
    P.nFreqPerBand = 9;

    % 动态 LTR 允许裕度。注意这不是总 LTR=1，而是留给动态晃动的部分。
    P.Mdyn = 0.60;

    % 仿真设置
    P.dt = 0.01;
    P.Tsim = 60;
    P.evalStart = 12;  % 跳过前段平滑进入与暂态
end

function [A,B,C,D,info] = build_model_B_scale(P)
    % 简化对象 B：两个耦合二阶模态。
    % x = [phi; phi_dot; theta; theta_dot]
    % phi   : 车辆/罐体侧倾代理角
    % theta : 液体主晃动代理角
    % input : ay
    % output: dynamic LTR proxy
    wr = 2*pi*P.f_roll;
    wl = 2*pi*P.f_liq;
    zr = P.zeta_roll;
    zl = P.zeta_liq;
    kc = P.k_couple;

    % 侧倾通道受 ay 和液体 theta 激励；液体通道受 ay 和 phi 加速度代理激励。
    A = [0, 1, 0, 0;
        -wr^2, -2*zr*wr, kc*wr^2, 0;
         0, 0, 0, 1;
         0.18*wl^2, 0, -wl^2, -2*zl*wl];

    B = [0;
         0.40;
         0;
         1.00];

    C = [P.gain_phi_to_ltr, 0.0, P.gain_theta_to_ltr, 0.0];
    D = P.gain_ay_to_ltr;

    info.description = 'coupled roll-slosh linear model B';
end

function Gabs = freq_response_mag(A,B,C,D,w)
    n = size(A,1);
    Gabs = zeros(size(w));
    I = eye(n);
    for k = 1:numel(w)
        G = C*((1i*w(k)*I - A)\B) + D;
        Gabs(k) = abs(G);
    end
end

function r = smooth_ramp(t,t0,t1)
    s = min(max((t-t0)/(t1-t0),0),1);
    r = 10*s.^3 - 15*s.^4 + 6*s.^5;
end

function [x,z] = simulate_lti(A,B,C,D,t,u)
    dt = t(2)-t(1);
    nx = size(A,1);
    x = zeros(numel(t),nx);
    z = zeros(numel(t),1);
    f = @(xx,uu) A*xx + B*uu;
    for k = 1:numel(t)-1
        uk = u(k);
        k1 = f(x(k,:).', uk);
        k2 = f(x(k,:).' + 0.5*dt*k1, uk);
        k3 = f(x(k,:).' + 0.5*dt*k2, uk);
        k4 = f(x(k,:).' + dt*k3, uk);
        x(k+1,:) = (x(k,:).' + dt/6*(k1+2*k2+2*k3+k4)).';
    end
    for k = 1:numel(t)
        z(k) = C*x(k,:).' + D*u(k);
    end
end

function E = band_power_projection(u,t,bandHz,nFreqPerBand,evalStart)
    idx = t >= evalStart;
    u = u(idx) - mean(u(idx));
    tt = t(idx);
    tt = tt - tt(1);
    N = numel(u);
    fList = linspace(bandHz(1), bandHz(2), nFreqPerBand);
    E = 0;
    for j = 1:numel(fList)
        c = cos(2*pi*fList(j)*tt);
        s = sin(2*pi*fList(j)*tt);
        ac = 2/N * sum(u .* c);
        as = 2/N * sum(u .* s);
        % 单个正弦分量的平均功率约为 A^2/2。
        E = E + 0.5*(ac^2 + as^2);
    end
end

function y = rms_local(x)
    y = sqrt(mean(x(:).^2));
end

function print_case_table(cases, Pcrit, P)
    fprintf('\n\n========== 测试用例结果表 ==========' );
    fprintf('\n%-4s %-34s %-10s %-10s %-10s %-10s %-8s', ...
        'No.','Case','Band/Pc','Total/Pc','PeakLTR','RMS_LTR','Safe');
    for i = 1:numel(cases)
        fprintf('\n%-4d %-34s %-10.2f %-10.2f %-10.3f %-10.3f %-8s', ...
            i, cases(i).short, cases(i).ratioBand, cases(i).ratioTotal, ...
            cases(i).peakAbsLTR, cases(i).rmsLTR, string(cases(i).safe));
    end
    fprintf('\n\n判读规则：Band/Pc <= 1 是 Qband 理论的保守安全判据；PeakLTR <= %.2f 是仿真观察到的动态安全判据。\n', P.Mdyn);
end

function plot_frequency_story(freqHz, Gabs, P, fWorst, GmaxBand, Pcrit, Acrit)
    figure('Name','Qband Step 1: identify dangerous band','Color','w','Position',[80 80 1180 480]);
    hold on;
    yl = [0, max(Gabs)*1.15];
    patch([P.bandHz(1) P.bandHz(2) P.bandHz(2) P.bandHz(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.92 0.92 0.92], 'EdgeColor','none', 'FaceAlpha',0.8);
    plot(freqHz, Gabs, 'LineWidth', 2.0);
    xline(P.f_liq, '--', '0.5 Hz 主晃动频率', 'LineWidth',1.2);
    plot(fWorst, GmaxBand, 'o', 'MarkerSize',8, 'LineWidth',1.8);
    grid on; box on;
    xlabel('输入频率 f (Hz)');
    ylabel('|G(j2\pi f)|  动态LTR/(m/s^2)');
    title('第1步：找出横向加速度最容易激发动态 LTR 的频带');
    legend('Qband 约束频带','ay \rightarrow dynamic LTR 频响','主晃动频率','带内最危险点','Location','best');
    text(1.05, 0.90*yl(2), sprintf('带内最坏增益 = %.3f\n保守功率阈值 P_{crit}=%.3f\n对应单频幅值 A_{crit}=%.3f m/s^2', ...
        GmaxBand, Pcrit, Acrit), 'BackgroundColor','w', 'EdgeColor',[0.75 0.75 0.75]);
end

function plot_case_summary(cases, P)
    n = numel(cases);
    peakLTR = [cases.peakAbsLTR];
    ratioBand = [cases.ratioBand];

    figure('Name','Qband Step 2: case summary','Color','w','Position',[100 100 1200 520]);
    subplot(1,2,1);
    bar(ratioBand); hold on;
    yline(1.0, '--', 'Qband 阈值', 'LineWidth',1.4);
    grid on; box on;
    set(gca,'XTick',1:n,'XTickLabel',{cases.short},'XTickLabelRotation',20);
    ylabel('危险频带输入功率 / P_{crit}');
    title('第2步：每个输入的“危险频带能量”是否超限');

    subplot(1,2,2);
    bar(peakLTR); hold on;
    yline(P.Mdyn, '--', '+动态LTR裕度', 'LineWidth',1.4);
    grid on; box on;
    set(gca,'XTick',1:n,'XTickLabel',{cases.short},'XTickLabelRotation',20);
    ylabel('峰值 |dynamic LTR|');
    title('第3步：仿真输出是否越过动态 LTR 裕度');
end

function plot_time_story(t, cases, P)
    showIdx = [1 2 3 4];
    figure('Name','Qband Step 3: time-domain meaning','Color','w','Position',[120 80 1280 760]);

    subplot(2,1,1); hold on;
    for i = showIdx
        plot(t, cases(i).ay, 'LineWidth',1.3);
    end
    grid on; box on;
    xlabel('t (s)'); ylabel('a_y input (m/s^2)');
    title('输入横向加速度：看起来都是振荡，但危险性取决于是否集中在主晃动频带');
    legend({cases(showIdx).short}, 'Location','bestoutside');
    xlim([0, min(max(t),40)]);

    subplot(2,1,2); hold on;
    for i = showIdx
        plot(t, cases(i).z, 'LineWidth',1.4);
    end
    yline(P.Mdyn, '--', '+动态LTR裕度', 'LineWidth',1.3);
    yline(-P.Mdyn, '--', '-动态LTR裕度', 'LineWidth',1.3);
    grid on; box on;
    xlabel('t (s)'); ylabel('dynamic LTR');
    title('输出动态 LTR：带内功率超限时，晃动响应更容易持续放大并触及安全裕度');
    legend([{cases(showIdx).short}, {'+Mdyn','-Mdyn'}], 'Location','bestoutside');
    xlim([0, min(max(t),40)]);
end

function plot_band_vs_total(cases, Pcrit, P)
    n = numel(cases);
    figure('Name','Qband Step 4: band power vs total power','Color','w','Position',[140 120 1000 620]);
    hold on;
    for i = 1:n
        if cases(i).safe
            mk = 'o';
        else
            mk = 'x';
        end
        plot(cases(i).totalPower/Pcrit, cases(i).bandPower/Pcrit, mk, 'MarkerSize',11, 'LineWidth',2.0);
        text(cases(i).totalPower/Pcrit + 0.03, cases(i).bandPower/Pcrit + 0.03, sprintf('Case %d',i));
    end
    xline(1.0, ':', '总功率=Pcrit');
    yline(1.0, '--', 'Qband功率=Pcrit');
    grid on; box on;
    xlabel('总输入平均功率 / P_{crit}');
    ylabel('危险频带输入功率 / P_{crit}');
    title('第4步：总功率不等于危险功率；Qband 约束的是会激发主晃动的那部分能量');
    subtitle(sprintf('目标频带 %.2f--%.2f Hz，主频 %.2f Hz', P.bandHz(1), P.bandHz(2), P.f_liq));
    legend('安全用例','Location','best');
end
