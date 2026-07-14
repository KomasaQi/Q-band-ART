function diagOut = run_qband_diagnostic(stopTime)
%RUN_QBAND_DIAGNOSTIC Run Q-band controller once and save compact diagnostics.
% This script does not save the Simulink model or full simOut.
    if nargin < 1
        stopTime = 60;
    end

    clear Pin5_MPC_Qband_TrajTrack_Ctrllor
    model = 'ALGO_TrajTrack_Qband';
    load_system([model '.slx']);
    blk = [model '/MPC_T1_SP_LTR_Cons_Ctrllor'];
    fprintf('OLD_FN=%s\n',get_param(blk,'FunctionName'));
    set_param(blk,'FunctionName','Pin5_MPC_Qband_TrajTrack_Ctrllor');
    fprintf('RUN_FN=%s\n',get_param(blk,'FunctionName'));
    set_param(model,'StopFcn','');
    set_param(model,'CloseFcn','');
    set_param(model,'Dirty','off');
    stopTimeStr = num2str(stopTime);
    set_param(model,'StopTime',stopTimeStr);

    simOut = sim(model,'StopTime',stopTimeStr);
    metricsQband = calc_5pin_segment_metrics(simOut,'qband');
    if exist('compare_target_segment_metrics.mat','file')
        S = load('compare_target_segment_metrics.mat');
        metricsTarget = S.metricsTarget;
    else
        metricsTarget = calc_5pin_segment_metrics('compare_target.mat','target');
        save('compare_target_segment_metrics.mat','metricsTarget');
    end

    compare = metricsQband(:,{'segment','peak_abs_LTR','max_lateral_error_m'});
    compare.target_peak_abs_LTR = metricsTarget.peak_abs_LTR;
    compare.target_max_lateral_error_m = metricsTarget.max_lateral_error_m;
    compare.ok_LTR = compare.peak_abs_LTR <= compare.target_peak_abs_LTR + 1e-9;
    compare.ok_err = compare.max_lateral_error_m <= compare.target_max_lateral_error_m + 1e-9;

    truck = simOut.TruckState;
    tTruck = truck.Time(:);
    truckData = truck.Data;
    x = truckData(:,9);
    y = truckData(:,10);
    ltr = interp_timeseries_local(simOut.LTR,tTruck);
    try
        ltrEq = interp_timeseries_local(simOut.LTR_eqval,tTruck);
    catch
        ltrEq = nan(size(ltr));
    end

    idxRoll = find(abs(ltr) >= 0.999,1,'first');
    if isempty(idxRoll)
        firstRollover = [nan nan nan nan];
    else
        firstRollover = [tTruck(idxRoll),x(idxRoll),y(idxRoll),ltr(idxRoll)];
    end

    try
        fallbackLog = evalin('base','qband_fallback_log');
    catch
        fallbackLog = zeros(0,5);
    end

    diagOut = struct();
    diagOut.metricsQband = metricsQband;
    diagOut.metricsTarget = metricsTarget;
    diagOut.compare = compare;
    diagOut.fallbackLog = fallbackLog;
    diagOut.fallbackCount = size(fallbackLog,1);
    diagOut.firstFallback = first_row_or_nan(fallbackLog,5);
    diagOut.lastFallback = last_row_or_nan(fallbackLog,5);
    diagOut.firstRollover = firstRollover; % [t, X, Y, LTR]
    diagOut.tEnd = tTruck(end);
    diagOut.xEnd = x(end);
    diagOut.yEnd = y(end);
    diagOut.finalTruckState = truckData(end,:);
    diagOut.finalLTR = ltr(end);
    diagOut.finalLTReq = ltrEq(end);
    diagOut.peakAbsLTR = max(abs(ltr),[],'omitnan');
    diagOut.peakAbsLTReq = max(abs(ltrEq),[],'omitnan');

    runStamp = datestr(now,'yyyymmdd_HHMMSS');
    runFile = ['qband_run_',runStamp,'.mat'];
    diagOut.runFile = runFile;
    save('qband_2dof_segment_metrics.mat','metricsQband','metricsTarget','compare');
    save(runFile,'simOut','metricsQband','metricsTarget','compare','diagOut','fallbackLog','-v7.3');
    save('qband_latest_run.mat','simOut','metricsQband','metricsTarget','compare','diagOut','fallbackLog','-v7.3');
    save('qband_diagnostic.mat','diagOut');
    disp('comparison:');
    disp(compare);
    fprintf('saved_run_file=%s\n',runFile);
    fprintf('saved_latest_file=qband_latest_run.mat\n');
    fprintf('fallback_count=%d\n',diagOut.fallbackCount);
    fprintf('first_fallback=[t X Y vx count]=%s\n',mat2str(diagOut.firstFallback,6));
    fprintf('last_fallback=[t X Y vx count]=%s\n',mat2str(diagOut.lastFallback,6));
    fprintf('first_rollover=[t X Y LTR]=%s\n',mat2str(diagOut.firstRollover,6));
    fprintf('end_state: t=%.3f, X=%.3f, Y=%.3f, final_LTR=%.4f, final_LTR_eq=%.4f\n', ...
        diagOut.tEnd,diagOut.xEnd,diagOut.yEnd,diagOut.finalLTR,diagOut.finalLTReq);
end

function y = interp_timeseries_local(ts,tTarget)
    t = ts.Time(:);
    data = ts.Data;
    data = data(:);
    [tUnique,ia] = unique(t,'stable');
    data = data(ia);
    y = interp1(tUnique,data,tTarget,'linear','extrap');
end

function row = first_row_or_nan(x,n)
    if isempty(x)
        row = nan(1,n);
    else
        row = x(1,:);
    end
end

function row = last_row_or_nan(x,n)
    if isempty(x)
        row = nan(1,n);
    else
        row = x(end,:);
    end
end
