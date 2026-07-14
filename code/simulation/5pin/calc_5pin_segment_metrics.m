function metrics = calc_5pin_segment_metrics(outOrFile, label)
%CALC_5PIN_SEGMENT_METRICS Segment metrics for the 5 independent DLC ranges.
%   Segments are defined by tractor X position:
%   1: 0-325 m, 2: 325-600 m, 3: 600-825 m,
%   4: 825-1025 m, 5: 1025-1200 m.

    if nargin < 2
        label = '';
    end
    if ischar(outOrFile) || isstring(outOrFile)
        S = load(outOrFile);
        if isfield(S,'out')
            out = S.out;
        elseif isfield(S,'simOut')
            out = S.simOut;
        else
            error('No variable named out or simOut found in %s.', outOrFile);
        end
    else
        out = outOrFile;
    end

    segEdges = [0 325; 325 600; 600 825; 825 1025; 1025 1200];
    nSeg = size(segEdges,1);
    [path,~,~] = getPath('new_5pin',1);
    pathX = path(:,1);
    pathY = path(:,2);

    truck = out.TruckState;
    tTruck = truck.Time(:);
    truckData = truck.Data;
    x = truckData(:,9);
    y = truckData(:,10);

    ltr = interp_timeseries(out.LTR,tTruck);
    if isprop(out,'LTR_eqval') || any(strcmp(out.who,'LTR_eqval'))
        ltrEq = interp_timeseries(out.LTR_eqval,tTruck);
    else
        ltrEq = nan(size(ltr));
    end

    eLat = local_path_error(x,y,pathX,pathY);

    metrics = table((1:nSeg)',segEdges(:,1),segEdges(:,2), ...
        nan(nSeg,1),nan(nSeg,1),nan(nSeg,1),nan(nSeg,1),nan(nSeg,1), ...
        'VariableNames',{'segment','x_start_m','x_end_m','peak_abs_LTR', ...
        'peak_abs_LTR_eqval','max_lateral_error_m','rms_lateral_error_m','sample_count'});

    for i = 1:nSeg
        if i < nSeg
            idx = x >= segEdges(i,1) & x < segEdges(i,2);
        else
            idx = x >= segEdges(i,1) & x <= segEdges(i,2);
        end
        idx = idx & isfinite(ltr) & isfinite(eLat);
        metrics.sample_count(i) = nnz(idx);
        if any(idx)
            metrics.peak_abs_LTR(i) = max(abs(ltr(idx)));
            metrics.peak_abs_LTR_eqval(i) = max(abs(ltrEq(idx)));
            metrics.max_lateral_error_m(i) = max(eLat(idx));
            metrics.rms_lateral_error_m(i) = rms(eLat(idx));
        end
    end

    if ~isempty(label)
        fprintf('\n%s segment metrics:\n',label);
    end
    disp(metrics);
end

function y = interp_timeseries(ts,tTarget)
    t = ts.Time(:);
    data = ts.Data;
    data = data(:);
    [tUnique,ia] = unique(t,'stable');
    data = data(ia);
    y = interp1(tUnique,data,tTarget,'linear','extrap');
end

function err = local_path_error(x,y,pathX,pathY)
    n = numel(x);
    err = nan(n,1);
    j0 = 1;
    for k = 1:n
        if ~isfinite(x(k)) || ~isfinite(y(k))
            continue;
        end
        lo = max(1,j0-500);
        hi = min(numel(pathX),j0+1200);
        [d2,jRel] = min((pathX(lo:hi)-x(k)).^2 + (pathY(lo:hi)-y(k)).^2);
        j0 = lo+jRel-1;
        err(k) = sqrt(d2);
    end
end
