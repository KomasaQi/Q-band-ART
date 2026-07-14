%% visualize_anti_exitation_5pins_adapt.m
% Regenerate the retained paper figure from saved ablation results.
% First run anti_exitation_nmpc_5pins_adapt.m once to create:
% anti_exitation_nmpc_5pins_adapt_results.mat

clear; clc; close all;

matlabDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(matlabDir, 'anti_exitation_nmpc_5pins_adapt_results.mat');
if ~exist(dataFile, 'file')
    error('Result file not found: %s\nRun anti_exitation_nmpc_5pins_adapt.m first.', dataFile);
end

S = load(dataFile, 'results', 'P');
results = S.results;
P = S.P;

[results, caseLabels] = remove_excluded_ablation_cases(results);
fprintf('Loaded %d ablation cases for visualization.\n', numel(results));

plot_segmentwise_statistics(results, P, caseLabels);
fprintf('Only the segment-wise paper figure is regenerated: Article-IEEE-ICVES-2026/fig/results3.png\n');

function [resultsOut, labels] = remove_excluded_ablation_cases(resultsIn)
    excluded = ["+ ay + dAy", "+ preview+feedback adaptive full"];
    keep = true(numel(resultsIn), 1);
    for i = 1:numel(resultsIn)
        keep(i) = ~any(string(resultsIn{i}.name) == excluded);
    end
    resultsOut = resultsIn(keep);

    labels = strings(numel(resultsOut), 1);
    for i = 1:numel(resultsOut)
        labels(i) = short_case_label(string(resultsOut{i}.name));
    end
end

function label = short_case_label(name)
    switch name
        case "Baseline: tracking only"
            label = "Base";
        case "+ ay hard constraint"
            label = "AyLim";
        case "+ ay cost only"
            label = "Ay";
        case "+ dAy only"
            label = "dAy";
        case "+ ddAy only"
            label = "ddAy";
        case "+ Qband only"
            label = "Qband";
        case "+ ay + dAy + ddAy"
            label = "Ay+dAy+ddAy";
        case "+ Qband soft only"
            label = "QSoft";
        case "+ fixed Qband full"
            label = "Fixed-Q";
        case "+ feedback adaptive Qband full"
            label = "Fb-Q";
        case "+ preview adaptive Qband full"
            label = "Pre-Q";
        case "+ fixed Qband + band soft"
            label = "Fixed-Q+Soft";
        case "+ preview adaptive + band soft"
            label = "Pre-Q+Soft";
        otherwise
            label = erase(name, "+ ");
    end
end

function plot_segmentwise_statistics(results, P, labels)
    if isempty(results) || isempty(results{1}.metrics.segmentRmsError)
        warning('No segment-wise metrics are available.');
        return;
    end

    rollover = getfield_with_default(P.ltr, 'rollover', 1.0);
    n = numel(results);
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
    rolloverFlags = segLTR >= rollover;

    fig = figure('Name','Segment-wise ablation statistics', ...
        'Color','w','Position',[90 70 1350 900]);
    tiledlayout(3, 1, 'TileSpacing','compact', 'Padding','compact');

    ax1 = nexttile;
    bar(segRms, 'LineWidth', 0.8);
    box on; grid on;
    ylabel('RMS error (m)');
    title('Tracking accuracy by independent DLC');
    set(ax1, 'XTick', 1:nSeg, 'XTickLabel', segNames);
    legend(labels, 'Location','northoutside', 'Orientation','horizontal', 'NumColumns', 5);

    ax2 = nexttile;
    bar(segPeak, 'LineWidth', 0.8);
    box on; grid on;
    ylabel('Peak |error| (m)');
    title('Peak tracking error by independent DLC');
    set(ax2, 'XTick', 1:nSeg, 'XTickLabel', segNames);

    ax3 = nexttile;
    bar(min(segLTR, rollover), 'LineWidth', 0.8);
    hold on; box on; grid on;
    [row, col] = find(rolloverFlags);
    if ~isempty(row)
        groupWidth = min(0.8, n/(n + 1.5));
        x = row - groupWidth/2 + (2*col-1) * groupWidth / (2*n);
        plot(x, rollover*ones(size(x)), 'rx', 'LineWidth', 1.4, 'MarkerSize', 7);
    end
    yline(rollover, 'k--', 'Rollover occurred', ...
        'LineWidth', 1.3, 'LabelHorizontalAlignment','left');
    ylim([0, 1.08*rollover]);
    ylabel('Peak |LTR_{est}|');
    xlabel('Independent DLC condition');
    title('Segment-wise rollover-risk proxy');
    set(ax3, 'XTick', 1:nSeg, 'XTickLabel', segNames);

    sgtitle('Segment-wise 5-pins validation statistics');
    save_paper_figure(fig, 'results3.png');
end

function value = getfield_with_default(S, fieldName, defaultValue)
    if isfield(S, fieldName)
        value = S.(fieldName);
    else
        value = defaultValue;
    end
end

function save_paper_figure(fig, fileName)
    scriptPath = mfilename('fullpath');
    matlabDir = fileparts(scriptPath);
    repoRoot = fileparts(fileparts(matlabDir));
    figDir = fullfile(repoRoot, 'Article-IEEE-ICVES-2026', 'fig');
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
    outPath = fullfile(figDir, fileName);
    try
        exportgraphics(fig, outPath, 'Resolution', 300);
    catch
        saveas(fig, outPath);
    end
    fprintf('Saved figure: %s\n', outPath);
end
