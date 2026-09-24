function compare_sleep_window_vs_MI_data()
% compare_sleep_window_vs_MI_data
% Compare mean sleep MI from sliding-window summary with MI_stat file.
%
% Writes a small MAT file with the numeric comparison.

try
    % Resolve repo root (assume this script sits in CoI-pipeline/Helpers)
    thisFile = mfilename('fullpath');
    repoRoot = fileparts(fileparts(thisFile));

    % Paths
    summaryMat = fullfile(repoRoot, 'Data', 'MI_Data', 'tmp_R060721', 'sleep_window_scan', 'R060721_BSLEEP_sleep_window_MI_summary.mat');
    dataMat = fullfile(repoRoot, 'Data', 'MI_Data', 'R060721sleepBSLEEP_MI_data.mat');

    if ~exist(summaryMat, 'file')
        error('Summary MAT not found: %s', summaryMat);
    end
    if ~exist(dataMat, 'file')
        error('Data MAT not found: %s', dataMat);
    end

    S = load(summaryMat, 'summary');
    if ~isfield(S, 'summary')
        error('Summary variable not found in %s', summaryMat);
    end
    summary = S.summary;

    % Collect mean(sleep_MI) for each window
    wr = summary.window_results;
    nW = numel(wr);
    window_means = nan(1, nW);
    for i = 1:nW
        si = wr(i).sleep_MI;
        window_means(i) = mean(si(:), 'omitnan');
    end
    slide_mean_mean = mean(window_means, 'omitnan');

    % Load MI_stat from the main MI data file
    D = load(dataMat, 'MI_stat');
    if ~isfield(D, 'MI_stat')
        error('MI_stat not found in %s', dataMat);
    end
    MI_stat = D.MI_stat;

    subj = 'R060721';
    cond = 'BSLEEP';
    chan_field = 'E14';

    if isfield(MI_stat, subj) && isfield(MI_stat.(subj), cond) && isfield(MI_stat.(subj).(cond).MI, chan_field)
        mi_arr = MI_stat.(subj).(cond).MI.(chan_field);
        file_mean = mean(mi_arr(:), 'omitnan');
    else
        error('Expected field MI_stat.%s.%s.MI.%s not found', subj, cond, chan_field);
    end

    diff = slide_mean_mean - file_mean;

    fprintf('\nComparison results:\n');
    fprintf(' Sliding-window mean of mean(sleep_MI) = %g\n', slide_mean_mean);
    fprintf(' MI_stat.%s.%s.MI.%s mean = %g\n', subj, cond, chan_field, file_mean);
    fprintf(' Difference (slide - file) = %g\n', diff);

    % Save result
    out.result.slide_mean_mean = slide_mean_mean;
    out.result.file_mean = file_mean;
    out.result.diff = diff;
    out.window_means = window_means;
    out_file = fullfile(repoRoot, 'CoI-pipeline', 'Helpers', 'compare_sleep_window_vs_MI_data_result.mat');
    save(out_file, '-struct', 'out');
    fprintf('Saved comparison to %s\n', out_file);

catch ME
    fprintf('compare_sleep_window_vs_MI_data failed: %s\n', ME.message);
    rethrow(ME);
end
end
