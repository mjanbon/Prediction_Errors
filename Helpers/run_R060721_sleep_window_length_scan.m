function run_R060721_sleep_window_length_scan(subject, condition_tag)
% run_R060721_sleep_window_length_scan
%
% Scan sleep-only MI across windows of increasing length for multiple
% channels.
% For each window length, this helper slides a window across the sleep
% trials, computes cnm_MI_stimtime on each window, and stores the mean MI
% across all windows of that length.
%
% Outputs:
%   - A MAT summary with the per-length results.
%   - A figure with 6 panels showing mean MI vs window length for E14,
%     E12, E10, E7, E4, and E1, each with the full-sleep mean as a dotted
%     reference line and the random-sampling mean as a red line.

try
    repoRoot = resolve_repo_root();
    addpath(genpath(repoRoot));

    topRoot = fileparts(repoRoot);
    cnmRoot = fullfile(topRoot, 'CNM');
    gcmiRoot = fullfile(topRoot, 'GCMI_master');
    if exist(cnmRoot, 'dir'), addpath(genpath(cnmRoot)); end
    if exist(gcmiRoot, 'dir'), addpath(genpath(gcmiRoot)); end

    rng(0, 'twister');
    USING_HPC = 2 %2 for queen mary, 0 for local

    if nargin < 1 || isempty(subject)
        subject = 'R060721';
    end
    if nargin < 2 || isempty(condition_tag)
        condition_tag = 'BSLEEP';
    end
    channel_list = [14, 12, 10, 7, 4, 1];

    % Window-length scan settings.
    min_window_len = [];
    max_window_len = [];
    window_len_step = [];
    max_windows_per_length = 25; % cap for speed; set [] to keep all windows.
    max_random_samples_per_length = 25; % cap for the random-sampling baseline.

    [basefold, datatype, ~, ~, ~, ~, ~, srate, activity_tag_1, ...
        ~, ~, activity_tag_2, ...
        deviant_group_number_2, standard_group_number_2, corrected, ...
        ~, baseline, ~, ~, ~] = Max_get_comparison_param(USING_HPC, 0);

    overVar_file = fullfile(basefold, datatype, sprintf('%s_%s.mat', subject, condition_tag));
    if ~exist(overVar_file, 'file')
        error('LFP file not found: %s', overVar_file);
    end

    fprintf('Using LFP file: %s\n', overVar_file);
    fprintf('Comparison tags: %s vs %s\n', activity_tag_1, activity_tag_2);

    [dvt_sleep, std_sleep] = load_trials_from_group_hyper(overVar_file, ...
        deviant_group_number_2, standard_group_number_2, corrected, srate);

    if any(channel_list > dvt_sleep.nbchan)
        error('One or more requested channels exceed available channels (%d).', dvt_sleep.nbchan);
    end

    % Sleep-only comparison baseline is computed exactly as in the existing helper.
    full_sleep_pair = baseline_normalize_pair(dvt_sleep, std_sleep, baseline);

    % Determine scan range.
    sleep_trial_count = std_sleep.trials;
    if isempty(min_window_len), min_window_len = max(50, floor(sleep_trial_count / 20)); end
    if isempty(max_window_len), max_window_len = sleep_trial_count; end
    if isempty(window_len_step), window_len_step = max(1, floor((max_window_len - min_window_len) / 10)); end

    window_lengths = min_window_len:window_len_step:max_window_len;
    if isempty(window_lengths) || window_lengths(end) ~= max_window_len
        window_lengths(end + 1) = max_window_len;
    end
    window_lengths = unique(window_lengths, 'stable');

    outdir = fullfile(repoRoot, 'Data', 'MI_Data', sprintf('tmp_%s', subject), 'sleep_window_length_scan');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    summary = struct();
    summary.subject = subject;
    summary.condition = condition_tag;
    summary.activity_tag_1 = activity_tag_1;
    summary.activity_tag_2 = activity_tag_2;
    summary.min_window_len = min_window_len;
    summary.max_window_len = max_window_len;
    summary.window_len_step = window_len_step;
    summary.window_lengths = window_lengths;
    summary.channel_list = channel_list;
    summary.max_random_samples_per_length = max_random_samples_per_length;

    fprintf('Scanning %d window lengths for %d channels...\n', numel(window_lengths), numel(channel_list));

    channel_results = repmat(struct( ...
        'channel', [], ...
        'full_sleep_MI', [], ...
        'full_sleep_sigMask', [], ...
        'full_sleep_mean', [], ...
        'results', []), 1, numel(channel_list));

    for iCh = 1:numel(channel_list)
        summary_channel = channel_list(iCh);
        chan_name = sprintf('E%d', summary_channel);

        sleep_data_full = extract_channel_trials(full_sleep_pair, summary_channel);
        stim_full = build_stim(size(sleep_data_full.dev, 1), size(sleep_data_full.std, 1));

        fprintf('Computing MI on full sleep data for %s...\n', chan_name);
        [full_sleep_MI, full_sleep_sig] = cnm_MI_stimtime([sleep_data_full.dev; sleep_data_full.std], stim_full, 0);
        full_sleep_mean = mean(full_sleep_MI(:), 'omitnan');

        fprintf('Scanning %d window lengths for %s...\n', numel(window_lengths), chan_name);

        results = repmat(struct( ...
            'window_len', [], ...
            'start_indices', [], ...
            'window_means', [], ...
            'slide_mean_mean', [], ...
            'random_indices', [], ...
            'random_means', [], ...
            'random_mean_mean', [], ...
            'n_windows', [], ...
            'n_random', [], ...
            'window_MI', []), 1, numel(window_lengths));

        for iLen = 1:numel(window_lengths)
            win_len = window_lengths(iLen);
            max_start = sleep_trial_count - win_len + 1;
            if max_start < 1
                error('Window length %d exceeds available sleep trials (%d).', win_len, sleep_trial_count);
            end

            if isempty(max_windows_per_length) || max_windows_per_length >= max_start
                start_indices = 1:max_start;
            else
                start_indices = unique(round(linspace(1, max_start, max_windows_per_length)));
            end

            nRandom = max_start;
            if ~isempty(max_random_samples_per_length)
                nRandom = min(nRandom, max_random_samples_per_length);
            end

            nStarts = numel(start_indices);
            window_means = nan(1, nStarts);
            window_MI = cell(1, nStarts);
            random_indices = cell(1, nRandom);
            random_means = nan(1, nRandom);

            fprintf('  %s length %d/%d: win_len=%d, windows=%d\n', chan_name, iLen, numel(window_lengths), win_len, nStarts);

            for w = 1:nStarts
                start_idx = start_indices(w);
                stop_idx = start_idx + win_len - 1;

                sleep_window = trim_window_pair(dvt_sleep, std_sleep, start_idx, stop_idx);
                sleep_window = baseline_normalize_pair(sleep_window.dev, sleep_window.std, baseline);
                sleep_data = extract_channel_trials(sleep_window, summary_channel);
                stim = build_stim(size(sleep_data.dev, 1), size(sleep_data.std, 1));

                [sleep_MI, ~] = cnm_MI_stimtime([sleep_data.dev; sleep_data.std], stim, 0);
                window_MI{w} = sleep_MI;
                window_means(w) = mean(sleep_MI(:), 'omitnan');
            end

            % Random-sampling baseline: sample win_len trials across the
            % entire sleep set rather than taking a contiguous block.
            fprintf('    %s length %d/%d: random samples=%d\n', chan_name, iLen, numel(window_lengths), nRandom);
            for r = 1:nRandom
                rand_idx = sort(randperm(sleep_trial_count, win_len));
                random_indices{r} = rand_idx;

                random_window = subset_window_pair(dvt_sleep, std_sleep, rand_idx);
                random_window = baseline_normalize_pair(random_window.dev, random_window.std, baseline);
                random_data = extract_channel_trials(random_window, summary_channel);
                stim = build_stim(size(random_data.dev, 1), size(random_data.std, 1));

                [random_MI, ~] = cnm_MI_stimtime([random_data.dev; random_data.std], stim, 0);
                random_means(r) = mean(random_MI(:), 'omitnan');
            end

            results(iLen).window_len = win_len;
            results(iLen).start_indices = start_indices;
            results(iLen).window_means = window_means;
            results(iLen).slide_mean_mean = mean(window_means, 'omitnan');
            results(iLen).random_indices = random_indices;
            results(iLen).random_means = random_means;
            results(iLen).random_mean_mean = mean(random_means, 'omitnan');
            results(iLen).n_windows = nStarts;
            results(iLen).n_random = nRandom;
            results(iLen).window_MI = window_MI;
        end

        channel_results(iCh).channel = summary_channel;
        channel_results(iCh).full_sleep_MI = full_sleep_MI;
        channel_results(iCh).full_sleep_sigMask = full_sleep_sig;
        channel_results(iCh).full_sleep_mean = full_sleep_mean;
        channel_results(iCh).results = results;
        channel_results(iCh).slide_mean_means = arrayfun(@(r) r.slide_mean_mean, results);
        channel_results(iCh).random_mean_means = arrayfun(@(r) r.random_mean_mean, results);
    end

    summary.channel_results = channel_results;
    summary.full_sleep_means = arrayfun(@(r) r.full_sleep_mean, channel_results);

    save(fullfile(outdir, sprintf('%s_%s_sleep_window_length_scan.mat', subject, condition_tag)), ...
        'summary', '-v7.3');

    % Plot: 6 panels, one per channel.
    figure('Position', [50 50 1500 900]);
    tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    for iCh = 1:numel(channel_list)
        chan_name = sprintf('E%d', channel_list(iCh));
        nexttile;
        plot(window_lengths, channel_results(iCh).slide_mean_means, '-o', 'LineWidth', 2, 'MarkerSize', 5); hold on;
        plot(window_lengths, channel_results(iCh).random_mean_means, '-o', 'Color', [0.85 0 0], 'LineWidth', 2, 'MarkerSize', 5);
        yline(channel_results(iCh).full_sleep_mean, 'k:', 'LineWidth', 2);
        xlabel('Sleep window length (trials)');
        ylabel('Mean MI across windows');
        title(sprintf('%s', chan_name));
        legend({'slide\_mean\_mean', 'random sample mean', 'full sleep mean'}, 'Location', 'best');
        grid on;
    end

    sgtitle(sprintf('%s %s sleep window-length scan', subject, condition_tag));

    saveas(gcf, fullfile(outdir, sprintf('%s_%s_sleep_window_length_scan.png', subject, condition_tag)));
    saveas(gcf, fullfile(outdir, sprintf('%s_%s_sleep_window_length_scan.fig', subject, condition_tag)));

    fprintf('\nSaved summary to %s\n', fullfile(outdir, sprintf('%s_%s_sleep_window_length_scan.mat', subject, condition_tag)));
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

function window_pair = subset_window_pair(dvt, std, trial_indices)
window_pair = struct();
window_pair.dev = dvt;
window_pair.std = std;
window_pair.dev.data = dvt.data(:,:,trial_indices);
window_pair.std.data = std.data(:,:,trial_indices);
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

function data = extract_channel_trials(window_pair, summary_channel)
tmp = squeeze(window_pair.dev.data(summary_channel, :, :));
data.dev = tmp';
tmp = squeeze(window_pair.std.data(summary_channel, :, :));
data.std = tmp';
end

function stim = build_stim(nDev, nStd)
stim = [zeros(1, nDev), ones(1, nStd)]';
end
