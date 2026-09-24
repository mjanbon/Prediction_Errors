function run_R060721_sleep_randomsample_test()
% Compare full sleep MI against first-N and random-N sleep subsets.
% This is a diagnostic helper, not a pipeline replacement.

try
    cwd = pwd;
    if endsWith(cwd, filesep)
        cwd = cwd(1:end-1);
    end
    if endsWith(cwd, 'CoI-pipeline')
        repoRoot = cwd;
    else
        repoRoot = fullfile(cwd, 'CoI-pipeline');
    end
    addpath(genpath(repoRoot));
    topRoot = fileparts(repoRoot);
    cnmRoot = fullfile(topRoot, 'CNM');
    gcmiRoot = fullfile(topRoot, 'GCMI_master');
    if exist(cnmRoot, 'dir'), addpath(genpath(cnmRoot)); end
    if exist(gcmiRoot, 'dir'), addpath(genpath(gcmiRoot)); end

    subject = 'R060721';
    condition_tag = 'BSLEEP';
    chan = 14;
    nRepeats = 20;
    rng(0);

    overVar_file = fullfile('C:','Users','maxim','OneDrive','Desktop','PhD','Juho Code','Co-I Pipeline','CoI-pipeline','Data','Drosophila_LFP', sprintf('%s_%s.mat', subject, condition_tag));
    if ~exist(overVar_file, 'file')
        error('Requested LFP file not found: %s', overVar_file);
    end

    [~, ~, ~, ~, ~, ~, ~, srate, ~, deviant_group_number, standard_group_number, corrected, ~, baseline, ~, ~, ~] = Max_get_param(0, 0);
    [~, ~, ~, ~, ~, ~, ~, srate2, ~, dev1, std1, ~, dev2, std2, corrected2, ~, baseline2, ~, ~, ~] = Max_get_comparison_param(0, 0);

    [dvt_sleep, std_sleep] = load_trials_from_group_hyper(overVar_file, deviant_group_number, standard_group_number, corrected, srate);
    [dvt_wake, std_wake] = load_trials_from_group_hyper(overVar_file, dev1, std1, corrected2, srate2);
    [dvt_sleep_cmp, std_sleep_cmp] = load_trials_from_group_hyper(overVar_file, dev2, std2, corrected2, srate2);

    [sleep_data_full, sleep_stim_full] = prep_condition(dvt_sleep, std_sleep, baseline, chan);
    [wake_data_full, wake_stim_full] = prep_condition(dvt_wake, std_wake, baseline2, chan);
    [sleep_cmp_data_full, sleep_cmp_stim_full] = prep_condition(dvt_sleep_cmp, std_sleep_cmp, baseline2, chan);

    nMatch = min(size(wake_data_full, 1), size(sleep_cmp_data_full, 1));

    % Keep class balance when building matched subsets: selecting 1:nMatch
    % would only grab the first class block and can force MI to zero.
    wake_class0 = sum(wake_stim_full == 0);
    wake_class1 = sum(wake_stim_full == 1);
    sleep_class0_idx = find(sleep_stim_full == 0);
    sleep_class1_idx = find(sleep_stim_full == 1);

    n0 = min(wake_class0, numel(sleep_class0_idx));
    n1 = min(wake_class1, numel(sleep_class1_idx));
    nMatch = n0 + n1;
    firstN_idx = [sleep_class0_idx(1:n0); sleep_class1_idx(1:n1)];

    fprintf('Full sleep trials: %d\n', size(sleep_data_full, 1));
    fprintf('Matched subset size: %d\n', nMatch);

    full_sleep_MI = cnm_MI_stimtime(sleep_data_full, sleep_stim_full, 0);
    firstN_sleep_MI = cnm_MI_stimtime(sleep_data_full(firstN_idx, :), sleep_stim_full(firstN_idx), 0);

    rmse_firstN = sqrt(mean((double(firstN_sleep_MI(:)) - double(full_sleep_MI(:))).^2, 'omitnan'));
    fprintf('full_sleep_MI: size=%s mean=%g max=%g\n', mat2str(size(full_sleep_MI)), mean(full_sleep_MI(:), 'omitnan'), max(full_sleep_MI(:)));
    fprintf('firstN_sleep_MI: size=%s mean=%g max=%g RMSE_vs_full=%g\n', mat2str(size(firstN_sleep_MI)), mean(firstN_sleep_MI(:), 'omitnan'), max(firstN_sleep_MI(:)), rmse_firstN);

    random_rmse = nan(nRepeats, 1);
    random_mean = nan(nRepeats, 1);
    best_rmse = inf;
    best_mi = [];
    best_idx = [];

    for r = 1:nRepeats
        idx0 = sleep_class0_idx(randperm(numel(sleep_class0_idx), n0));
        idx1 = sleep_class1_idx(randperm(numel(sleep_class1_idx), n1));
        idx = [idx0; idx1];
        rand_mi = cnm_MI_stimtime(sleep_data_full(idx, :), sleep_stim_full(idx), 0);
        random_mean(r) = mean(rand_mi(:), 'omitnan');
        random_rmse(r) = sqrt(mean((double(rand_mi(:)) - double(full_sleep_MI(:))).^2, 'omitnan'));
        if random_rmse(r) < best_rmse
            best_rmse = random_rmse(r);
            best_mi = rand_mi;
            best_idx = idx;
        end
        fprintf('repeat %02d/%02d: mean=%g RMSE_vs_full=%g\n', r, nRepeats, random_mean(r), random_rmse(r));
    end

    fprintf('\n--- Summary ---\n');
    fprintf('full_sleep_MI mean=%g max=%g\n', mean(full_sleep_MI(:), 'omitnan'), max(full_sleep_MI(:)));
    fprintf('firstN_sleep_MI mean=%g max=%g RMSE_vs_full=%g\n', mean(firstN_sleep_MI(:), 'omitnan'), max(firstN_sleep_MI(:)), rmse_firstN);
    fprintf('random subset RMSE: mean=%g min=%g max=%g\n', mean(random_rmse, 'omitnan'), min(random_rmse), max(random_rmse));
    fprintf('best random subset RMSE_vs_full=%g\n', best_rmse);
    if ~isempty(best_idx)
        fprintf('best random subset first five indices: %s\n', mat2str(best_idx(1:min(5, numel(best_idx)))));
    end
    if ~isempty(best_mi)
        fprintf('best_random_sleep_MI mean=%g max=%g\n', mean(best_mi(:), 'omitnan'), max(best_mi(:)));
    end

    save(fullfile(repoRoot, 'Data', 'MI_Data', 'tmp_R060721', 'R060721_sleep_randomsample_test.mat'), ...
        'full_sleep_MI', 'firstN_sleep_MI', 'random_rmse', 'random_mean', 'best_mi', 'best_idx', 'nMatch');
    fprintf('Saved results to %s\n', fullfile(repoRoot, 'Data', 'MI_Data', 'tmp_R060721', 'R060721_sleep_randomsample_test.mat'));
    fprintf('Done.\n');
catch ME
    fprintf('Script failed: %s\n', ME.message);
    for k = 1:length(ME.stack)
        s = ME.stack(k);
        fprintf('  %s:%d\n', s.file, s.line);
    end
    rethrow(ME);
end
end

function [data_mat, stim] = prep_condition(dvt, std, baseline, chan)
    data_dvt = permute(dvt.data, [1 3 2]);
    mean_data_dvt = squeeze(mean(data_dvt(:,:,baseline), 3));
    dvt.data = data_dvt - repmat(mean_data_dvt, [1 1 size(dvt.data, 2)]);
    bb1_dev_bl = permute(dvt.data, [1 3 2]);

    data_std = permute(std.data, [1 3 2]);
    mean_data_std = squeeze(mean(data_std(:,:,baseline), 3));
    std.data = data_std - repmat(mean_data_std, [1 1 size(std.data, 2)]);
    bb1_std_bl = permute(std.data, [1 3 2]);

    bb_dev = permute(bb1_dev_bl, [3 1 2]);
    bb_std = permute(bb1_std_bl, [3 1 2]);

    data_mat = [squeeze(bb_dev(:, chan, :)); squeeze(bb_std(:, chan, :))];
    stim = [zeros(size(bb_dev, 1), 1); ones(size(bb_std, 1), 1)];
end