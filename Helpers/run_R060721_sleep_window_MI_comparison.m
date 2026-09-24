function run_R060721_sleep_window_MI_comparison()
% run_R060721_sleep_window_MI_comparison
%
% Sliding-window helper for the MI comparison pipeline.
% It keeps the wake trials fixed at the first N trials, then slides a window
% of length N across the sleep trials: sleep window = n : n+N-1.
%
% This mirrors the trial-cutting logic used in the comparison pipeline,
% but replaces the deterministic "first N" sleep slice with a moving window.
%
% Output:
%   Saves a .mat summary and a quick-look figure for the chosen channel.

try
    repoRoot = resolve_repo_root();
    addpath(genpath(repoRoot));

    topRoot = fileparts(repoRoot);
    cnmRoot = fullfile(topRoot, 'CNM');
    gcmiRoot = fullfile(topRoot, 'GCMI_master');
    if exist(cnmRoot, 'dir'), addpath(genpath(cnmRoot)); end
    if exist(gcmiRoot, 'dir'), addpath(genpath(gcmiRoot)); end

    rng(0, 'twister');

    subject = 'R060721';
    condition_tag = 'BSLEEP';
    summary_channel = 14; % E14 by default; change if you want another channel.
    threshold = 0.05;
    exploratory_kperm = 10
    ; % smaller than the full pipeline to keep the helper practical
    sleep_window_step = [];

    [basefold, datatype, all_con, ~, ~, ~, ~, srate, activity_tag_1, ...
        deviant_group_number_1, standard_group_number_1, activity_tag_2, ...
        deviant_group_number_2, standard_group_number_2, corrected, ...
        ~, baseline, ~, ~, kperm] = Max_get_comparison_param(0, 0);

    if isempty(exploratory_kperm)
        exploratory_kperm = kperm;
    end

    overVar_file = fullfile(basefold, datatype, sprintf('%s_%s.mat', subject, condition_tag));
    if ~exist(overVar_file, 'file')
        error('LFP file not found: %s', overVar_file);
    end

    fprintf('Using LFP file: %s\n', overVar_file);
    fprintf('Comparison tags: %s vs %s\n', activity_tag_1, activity_tag_2);

    % Load wake and sleep trial structs using the same group loading logic as the pipeline.
    [dvt_wake, std_wake] = load_trials_from_group_hyper(overVar_file, ...
        deviant_group_number_1, standard_group_number_1, corrected, srate);
    [dvt_sleep, std_sleep] = load_trials_from_group_hyper(overVar_file, ...
        deviant_group_number_2, standard_group_number_2, corrected, srate);

    % The original comparison code trims both conditions to the shorter std count.
    trialnum = min(std_wake.trials, std_sleep.trials);
    if trialnum < 1
        error('Not enough trials to build a comparison window.');
    end

    % Fixed wake slice: first N trials, matching the current pipeline behavior.
    wake_start = 1;
    wake_stop = trialnum;

    % Sleep windows slide across the longer sleep condition.
    max_sleep_start = std_sleep.trials - trialnum + 1;
    if max_sleep_start < 1
        error('Sleep data is shorter than the wake match length.');
    end

    if isempty(sleep_window_step)
        sleep_window_step = max(1, floor(trialnum / 10));
    end
    sleep_starts = 1:sleep_window_step:max_sleep_start;
    if sleep_starts(end) ~= max_sleep_start
        sleep_starts(end + 1) = max_sleep_start;
    end

    fprintf('Matched trial length N = %d\n', trialnum);
    fprintf('Sleep windows: %d starts, step=%d\n', length(sleep_starts), sleep_window_step);

    outdir = fullfile(repoRoot, 'Data', 'MI_Data', 'tmp_R060721', 'sleep_window_scan');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    if summary_channel > dvt_wake.nbchan
        error('Requested summary channel E%d exceeds available channels (%d).', summary_channel, dvt_wake.nbchan);
    end

    summary = struct();
    summary.subject = subject;
    summary.condition = condition_tag;
    summary.activity_tag_1 = activity_tag_1;
    summary.activity_tag_2 = activity_tag_2;
    summary.trialnum = trialnum;
    summary.sleep_starts = sleep_starts;
    summary.sleep_window_step = sleep_window_step;
    summary.summary_channel = summary_channel;
    summary.exploratory_kperm = exploratory_kperm;
    summary.threshold = threshold;

    wake_window = trim_window_pair(dvt_wake, std_wake, wake_start, wake_stop);
    wake_window = baseline_normalize_pair(wake_window.dev, wake_window.std, baseline);

    % Prepare fixed wake data (trials x samples) for full-sleep comparison
    tmp = squeeze(wake_window.dev.data(summary_channel, :, :));
    wake_dev_fixed = tmp';
    tmp = squeeze(wake_window.std.data(summary_channel, :, :));
    wake_std_fixed = tmp';
    wake_data_fixed = [wake_dev_fixed; wake_std_fixed];

    % Compute full-sleep MI (baseline-normalize using the same baseline)
    full_sleep_pair = baseline_normalize_pair(dvt_sleep, std_sleep, baseline);
    tmp = squeeze(full_sleep_pair.dev.data(summary_channel, :, :));
    full_dev = tmp';
    tmp = squeeze(full_sleep_pair.std.data(summary_channel, :, :));
    full_std = tmp';
    sleep_data_full = [full_dev; full_std];
    stim_full = [zeros(1, size(full_dev, 1)), ones(1, size(full_std, 1))]';

    fprintf('Computing MI on full sleep data for comparison (this may take a moment)...\n');
    [sleep_MI_full, sleep_sig_full] = cnm_MI_stimtime(sleep_data_full, stim_full, 0);

    % Compute MI for fixed wake using its own stim vector, then compute
    % simple difference to full-sleep MI. We skip the permutation test here
    % because `cnm_MI_stimtime_perm_diff` requires the same number of
    % trials in both inputs (matched stim length).
    stim_wake = [zeros(1, size(wake_dev_fixed, 1)), ones(1, size(wake_std_fixed, 1))]';
    [wake_MI_fixed, wake_sig_fixed] = cnm_MI_stimtime(wake_data_fixed, stim_wake, 0);
    MI_1_full = wake_MI_fixed;
    MI_2_full = sleep_MI_full;
    MI_diff_full = MI_2_full - MI_1_full;
    sigMask_full = nan(size(MI_diff_full));

    chan_name = sprintf('E%d', summary_channel);
    window_results = repmat(struct( ...
        'sleep_start', [], ...
        'sleep_stop', [], ...
        'cmp_MI_diff', [], ...
        'cmp_MI_1', [], ...
        'cmp_MI_2', [], ...
        'cmp_sigMask', [], ...
        'sleep_MI', [], ...
        'sleep_sigMask', []), 1, length(sleep_starts));

    fprintf('\nRunning sliding windows...\n');
    for w = 1:length(sleep_starts)
        sleep_start = sleep_starts(w);
        sleep_stop = sleep_start + trialnum - 1;

        sleep_window = trim_window_pair(dvt_sleep, std_sleep, sleep_start, sleep_stop);
        sleep_window = baseline_normalize_pair(sleep_window.dev, sleep_window.std, baseline);

        % Extract channel-specific data: original shape is [chan x samples x trials]
        % we want [trials x samples] for the MI functions.
        tmp = squeeze(wake_window.dev.data(summary_channel, :, :));
        wake_dev = tmp'; % trials x samples
        tmp = squeeze(wake_window.std.data(summary_channel, :, :));
        wake_std = tmp';
        tmp = squeeze(sleep_window.dev.data(summary_channel, :, :));
        sleep_dev = tmp';
        tmp = squeeze(sleep_window.std.data(summary_channel, :, :));
        sleep_std = tmp';

        wake_data = [wake_dev; wake_std];
        sleep_data = [sleep_dev; sleep_std];
        % stim labels: length must match number of trials (rows)
        stim = [zeros(1, size(wake_dev, 1)), ones(1, size(wake_std, 1))]';

        [MI_diff, MI_1, MI_2, sigMask] = cnm_MI_stimtime_perm_diff( ...
            wake_data, sleep_data, stim, exploratory_kperm, threshold);

        [sleep_MI, sleep_sig] = cnm_MI_stimtime(sleep_data, stim, 0);

        window_results(w).sleep_start = sleep_start;
        window_results(w).sleep_stop = sleep_stop;
        window_results(w).cmp_MI_diff = MI_diff;
        window_results(w).cmp_MI_1 = MI_1;
        window_results(w).cmp_MI_2 = MI_2;
        window_results(w).cmp_sigMask = sigMask;
        window_results(w).sleep_MI = sleep_MI;
        window_results(w).sleep_sigMask = sleep_sig;

        fprintf('  Window %d/%d: sleep trials %d:%d\n', w, length(sleep_starts), sleep_start, sleep_stop);
    end

    % Save a compact summary plus the full window results for later inspection.
    summary.window_results = window_results;
    summary.window_centers = sleep_starts + floor((trialnum - 1) / 2);

    % Attach full-sleep MI and comparison results for easy reference
    summary.full_sleep_MI = sleep_MI_full;
    summary.full_sleep_sigMask = sleep_sig_full;
    summary.full_sleep_MI_mean = mean(sleep_MI_full(:), 'omitnan');
    summary.full_cmp_MI_diff = MI_diff_full;
    summary.full_cmp_MI_1 = MI_1_full;
    summary.full_cmp_MI_2 = MI_2_full;
    summary.full_cmp_sigMask = sigMask_full;

    % If a saved MI_stat file exists, load and compare the saved MI for this
    % subject/condition/channel and store the difference in the summary.
    try
        miStatFile = fullfile(repoRoot, 'Data', 'MI_Data', sprintf('%ssleep%s_MI_data.mat', subject, condition_tag));
        if exist(miStatFile, 'file')
            D = load(miStatFile, 'MI_stat');
            if isfield(D, 'MI_stat')
                MI_stat = D.MI_stat;
                chan_field = chan_name; % e.g. 'E14'
                if isfield(MI_stat, subject) && isfield(MI_stat.(subject), condition_tag) && isfield(MI_stat.(subject).(condition_tag).MI, chan_field)
                    saved_mi = MI_stat.(subject).(condition_tag).MI.(chan_field);
                    saved_mi_mean = mean(saved_mi(:), 'omitnan');
                    summary.saved_MI_stat_mean = saved_mi_mean;
                    summary.full_vs_saved_mean_diff = summary.full_sleep_MI_mean - saved_mi_mean;
                    fprintf('\nFound MI_stat file: %s\n', miStatFile);
                    fprintf(' summary.full_sleep_MI_mean = %g\n', summary.full_sleep_MI_mean);
                    fprintf(' MI_stat.%s.%s.MI.%s mean = %g\n', subject, condition_tag, chan_field, saved_mi_mean);
                    fprintf(' Difference (full - saved) = %g\n', summary.full_vs_saved_mean_diff);
                else
                    fprintf('\nMI_stat does not contain expected field for %s/%s/%s.\n', subject, condition_tag, chan_field);
                end
            else
                fprintf('\nMI_stat variable not found in %s\n', miStatFile);
            end
        else
            fprintf('\nMI_stat file not found: %s\n', miStatFile);
        end
    catch ME2
        fprintf('\nWarning: failed to load/compare MI_stat file: %s\n', ME2.message);
    end

    save(fullfile(outdir, sprintf('%s_%s_sleep_window_MI_summary.mat', subject, condition_tag)), ...
        'summary', '-v7.3');

    % Create a quick-look plot for the summary channel using the mean MI over time.
    cmp_means = nan(1, length(window_results));
    sleep_means = nan(1, length(window_results));
    for w = 1:length(window_results)
        cmp_means(w) = mean(window_results(w).cmp_MI_diff(:), 'omitnan');
        sleep_means(w) = mean(window_results(w).sleep_MI(:), 'omitnan');
    end

    figure('Position', [100 100 1100 450]);
    plot(summary.window_centers, cmp_means, '-o', 'LineWidth', 2, 'DisplayName', 'comparison MI diff'); hold on;
    plot(summary.window_centers, sleep_means, '-o', 'LineWidth', 2, 'DisplayName', 'sleep-only MI');
    yline(0, '--k');
    xlabel('Sleep window center trial');
    ylabel('Mean MI (bits)');
    title(sprintf('%s %s sliding sleep windows, %s', subject, condition_tag, chan_name));
    legend('Location', 'best');
    grid on;

    saveas(gcf, fullfile(outdir, sprintf('%s_%s_sleep_window_MI_summary_%s.png', subject, condition_tag, chan_name)));
    saveas(gcf, fullfile(outdir, sprintf('%s_%s_sleep_window_MI_summary_%s.fig', subject, condition_tag, chan_name)));

    fprintf('\nSaved summary to %s\n', fullfile(outdir, sprintf('%s_%s_sleep_window_MI_summary.mat', subject, condition_tag)));
    fprintf('Done.\n');

catch ME
    fprintf('Helper failed: %s\n', ME.message);
    rethrow(ME);
end
end

function repoRoot = resolve_repo_root()
cwd = pwd;
if endsWith(cwd, filesep)
    cwd = cwd(1:end-1);
end
if endsWith(cwd, 'CoI-pipeline')
    repoRoot = cwd;
else
    repoRoot = fullfile(cwd, 'CoI-pipeline');
end
end

function window_pair = trim_window_pair(dvt, std, start_idx, stop_idx)
window_pair = struct();
window_pair.dev = dvt;
window_pair.std = std;
window_pair.dev.data = dvt.data(:,:,start_idx:stop_idx);
window_pair.std.data = std.data(:,:,start_idx:stop_idx);
window_pair.dev.trials = size(window_pair.dev.data, 3);
window_pair.std.trials = size(window_pair.std.data, 3);
end

function window_pair = baseline_normalize_pair(dev_struct, std_struct, baseline)
window_pair = struct('dev', dev_struct, 'std', std_struct);

% Deviant trials
data_dev = permute(window_pair.dev.data, [1 3 2]);
mean_dev = squeeze(mean(data_dev(:,:,baseline), 3));
window_pair.dev.data = data_dev - repmat(mean_dev, [1 1 size(window_pair.dev.data, 2)]);
window_pair.dev.data = permute(window_pair.dev.data, [1 3 2]);

% Standard trials
data_std = permute(window_pair.std.data, [1 3 2]);
mean_std = squeeze(mean(data_std(:,:,baseline), 3));
window_pair.std.data = data_std - repmat(mean_std, [1 1 size(window_pair.std.data, 2)]);
window_pair.std.data = permute(window_pair.std.data, [1 3 2]);
end
