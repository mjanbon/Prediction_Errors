function run_R060721_compare_pretrim_vs_single_E14()
% Run pretrim comparison pipeline and compare E14 MI_diff to
% single-condition MI(sleep) - MI(wake) using main_MI_ERP-style preprocessing.

try
    repoRoot = find_repo_root(pwd);

    addpath(genpath(repoRoot));
    topRoot = fileparts(repoRoot);
    cnmRoot = fullfile(topRoot, 'CNM');
    gcmiRoot = fullfile(topRoot, 'GCMI_master');
    if exist(cnmRoot, 'dir'), addpath(genpath(cnmRoot)); end
    if exist(gcmiRoot, 'dir'), addpath(genpath(gcmiRoot)); end

    subject = 'R060721';
    conditionTag = 'BSLEEP';
    chan = 14;

    fprintf('Repo root: %s\n', repoRoot);

    %% 1) Run pretrim comparison pipeline
    fprintf('\nRunning main_MI_ERP_comparison_pretrim...\n');
    main_MI_ERP_comparison_pretrim(1);

    %% 2) Load MI_diff(E14) from pretrim output
    miDataDir = fullfile(repoRoot, 'Data', 'MI_Data');
    pretrimFile = fullfile(miDataDir, sprintf('%s_pretrim_wake_sleep_%s_MI_data.mat', subject, conditionTag));
    if ~exist(pretrimFile, 'file')
        error('Pretrim MI file not found: %s', pretrimFile);
    end

    S = load(pretrimFile);
    if ~isfield(S, 'MI_stat') || ~isfield(S.MI_stat, 'MI_diff') || ~isfield(S.MI_stat.MI_diff, 'E14')
        error('Unable to find MI_stat.MI_diff.E14 in %s', pretrimFile);
    end
    pretrim_diff_E14 = double(S.MI_stat.MI_diff.E14(:));

    %% 3) Recompute single-condition sleep and wake MI (main_MI_ERP-style)
    fprintf('Recomputing sleep and wake MI with main_MI_ERP preprocessing (E14)...\n');

    [~, ~, ~, ~, ~, ~, ~, srateSleep, ~, devSleep, stdSleep, correctedSleep, ~, baselineSleep, ~, ~, ~] = Max_get_param(0, 0);
    [~, ~, ~, ~, ~, ~, ~, srateCmp, ~, devWake, stdWake, ~, devSleepCmp, stdSleepCmp, correctedCmp, ~, baselineCmp, ~, ~, ~] = Max_get_comparison_param(0, 0);

    lfpFile = fullfile(repoRoot, 'Data', 'Drosophila_LFP', sprintf('%s_%s.mat', subject, conditionTag));
    if ~exist(lfpFile, 'file')
        error('LFP file not found: %s', lfpFile);
    end

    [dvtSleep, stdSleepStruct] = load_trials_from_group_hyper(lfpFile, devSleep, stdSleep, correctedSleep, srateSleep);
    [dvtWake, stdWakeStruct] = load_trials_from_group_hyper(lfpFile, devWake, stdWake, correctedCmp, srateCmp);

    % Also load sleep groups from comparison params so the setup is explicit.
    [dvtSleepCmp, stdSleepCmpStruct] = load_trials_from_group_hyper(lfpFile, devSleepCmp, stdSleepCmp, correctedCmp, srateCmp);
    %#ok<NASGU> % Loaded for traceability; not directly used in the single-condition branch.

    sleep_E14 = compute_single_condition_mi(dvtSleep, stdSleepStruct, baselineSleep, chan);
    wake_E14 = compute_single_condition_mi(dvtWake, stdWakeStruct, baselineCmp, chan);

    single_diff_E14 = sleep_E14 - wake_E14;

    %% 4) Compare vectors
    if ~isequal(size(pretrim_diff_E14), size(single_diff_E14))
        error('Shape mismatch: pretrim=%s single=%s', mat2str(size(pretrim_diff_E14)), mat2str(size(single_diff_E14)));
    end

    delta = pretrim_diff_E14 - single_diff_E14;
    rmse = sqrt(mean(delta.^2, 'omitnan'));
    cc = corr(pretrim_diff_E14, single_diff_E14, 'rows', 'complete');

    fprintf('\n--- E14 Comparison ---\n');
    fprintf('pretrim MI_diff: mean=%g max=%g\n', mean(pretrim_diff_E14, 'omitnan'), max(pretrim_diff_E14));
    fprintf('single sleep-wake diff: mean=%g max=%g\n', mean(single_diff_E14, 'omitnan'), max(single_diff_E14));
    fprintf('delta (pretrim-single): mean=%g maxAbs=%g RMSE=%g corr=%g\n', mean(delta, 'omitnan'), max(abs(delta)), rmse, cc);

    outDir = fullfile(repoRoot, 'Data', 'MI_Data', 'tmp_R060721');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    outFile = fullfile(outDir, 'R060721_compare_pretrim_vs_single_E14.mat');
    save(outFile, 'pretrim_diff_E14', 'sleep_E14', 'wake_E14', 'single_diff_E14', 'delta', 'rmse', 'cc', 'pretrimFile');
    fprintf('Saved comparison output to %s\n', outFile);
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

function miVec = compute_single_condition_mi(dvt, stdStruct, baseline, chan)
% Reproduce main_MI_ERP preprocessing for one condition on one channel.

if dvt.trials < stdStruct.trials
    trialnum = dvt.trials;
    stdStruct.data = stdStruct.data(:, :, 1:trialnum);
    stdStruct.trials = trialnum;
else
    trialnum = stdStruct.trials;
    dvt.data = dvt.data(:, :, 1:trialnum);
    dvt.trials = trialnum;
end

data_dvt = permute(dvt.data, [1 3 2]);
mean_data_dvt = squeeze(mean(data_dvt(:, :, baseline), 3));
dvt.data = data_dvt - repmat(mean_data_dvt, [1 1 size(dvt.data, 2)]);
bb1_dev_bl = permute(dvt.data, [1 3 2]);

data_std = permute(stdStruct.data, [1 3 2]);
mean_data_std = squeeze(mean(data_std(:, :, baseline), 3));
stdStruct.data = data_std - repmat(mean_data_std, [1 1 size(stdStruct.data, 2)]);
bb1_std_bl = permute(stdStruct.data, [1 3 2]);

bb_dev = permute(bb1_dev_bl, [3 1 2]);
bb_std = permute(bb1_std_bl, [3 1 2]);

dataMat = [squeeze(bb_dev(:, chan, :)); squeeze(bb_std(:, chan, :))];
stim = [zeros(size(bb_dev, 1), 1); ones(size(bb_std, 1), 1)];

miVec = cnm_MI_stimtime(dataMat, stim, 0);
miVec = double(miVec(:));
end

function repoRoot = find_repo_root(startPath)
% Walk upward until a folder containing Code/ and Data/ is found.

repoRoot = startPath;
while true
    hasCode = exist(fullfile(repoRoot, 'Code'), 'dir') == 7;
    hasData = exist(fullfile(repoRoot, 'Data'), 'dir') == 7;
    if hasCode && hasData
        return;
    end

    parent = fileparts(repoRoot);
    if strcmp(parent, repoRoot)
        error('Could not locate CoI-pipeline root from start path: %s', startPath);
    end
    repoRoot = parent;
end
end
