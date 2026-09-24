function run_R060721_MI_comparison()
% run_R060721_MI_comparison.m
% Recompute single-condition MI (wake, sleep), comparison MI (wake vs sleep)
% for fly R060721, save intermediate inputs, and compare to legacy CoI outputs.
% Usage:
%   In MATLAB: run_R060721_MI_comparison
%   From shell: matlab -batch "run_R060721_MI_comparison"

try
    % Determine repo root robustly. If current working directory already
    % ends with 'CoI-pipeline', use it directly; otherwise append the folder.
    cwd = pwd;
    if endsWith(cwd, filesep) 
        cwd = cwd(1:end-1);
    end
    if endsWith(cwd, 'CoI-pipeline')
        repoRoot = cwd;
    else
        repoRoot = fullfile(cwd,'CoI-pipeline');
    end
    addpath(genpath(repoRoot));
    topRoot = fileparts(repoRoot);
    cnmRoot = fullfile(topRoot, 'CNM');
    gcmiRoot = fullfile(topRoot, 'GCMI_master');
    if exist(cnmRoot, 'dir'), addpath(genpath(cnmRoot)); end
    if exist(gcmiRoot, 'dir'), addpath(genpath(gcmiRoot)); end

    subject = 'R060721';
    condition_tag = 'BSLEEP';
    chan = 14; % electrode index for E14

    fprintf('Repo root: %s\n', repoRoot);

    %% Paths
    % Use the exact LFP path for R060721 as requested
    overVar_file = fullfile('C:','Users','maxim','OneDrive','Desktop','PhD','Juho Code','Co-I Pipeline','CoI-pipeline','Data','Drosophila_LFP', sprintf('%s_%s.mat', subject, condition_tag));
    if ~exist(overVar_file,'file')
        error('Requested LFP file not found: %s', overVar_file);
    end
    fprintf('Using LFP file: %s\n', overVar_file);
    outdir = fullfile(repoRoot, 'Data', 'MI_Data', 'tmp_R060721');
    if ~exist(outdir,'dir'), mkdir(outdir); end

    %% Get params for single-run pipeline
    fprintf('\nLoading single-run parameters...\n');
    try
        [basefold, datatype, all_con, condition, subject_list, participants, EoI, srate, activity_tag, deviant_group_number, standard_group_number, corrected, stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_param(0,0);
    catch
        warning('Max_get_param failed; falling back to default group indices');
        deviant_group_number = [3 4];
        standard_group_number = [9 10];
        corrected = 0; srate = 1000;
    end

    fprintf('Loading trials for single-run (sleep and wake) from %s\n', overVar_file);
    [dvt_s, std_s] = load_trials_from_group_hyper(overVar_file, deviant_group_number, standard_group_number, corrected, srate);

    % For wake groups, try standard comparison groups if available
    try
        [basefold2, datatype2, all_con2, condition2, subject2, participants2, EoI2, srate2, activity_tag_1, dev1, std1, activity_tag_2, dev2, std2, corrected2, stim_onset2, baseline2, start_cut_off2, end_cut_off2, kperm2] = Max_get_comparison_param(0,0);
    catch
        warning('Max_get_comparison_param failed; falling back to comparison group indices');
        dev1 = [7 8]; std1 = [11 12]; dev2 = [3 4]; std2 = [9 10]; corrected2 = corrected; srate2 = srate;
    end

    [dvt_w, std_w] = load_trials_from_group_hyper(overVar_file, dev1, std1, corrected2, srate2);
    [dvt_sl, std_sl] = load_trials_from_group_hyper(overVar_file, dev2, std2, corrected2, srate2);

    fprintf('Single-run trials: sleep dvt=%d std=%d; wake dvt=%d std=%d\n', dvt_s.trials, std_s.trials, dvt_w.trials, std_w.trials);
    fprintf('Comparison pipeline trials (before trimming): wake dvt=%d std=%d | sleep dvt=%d std=%d\n', dvt_w.trials, std_w.trials, dvt_sl.trials, std_sl.trials);

    %% Save copnormed data / stim used by downstream estimator for reproducibility
    try
        % Many pipeline functions expect 'dvt' and 'std' structs which include copnormed data
        save(fullfile(outdir, 'R060721_trial_structs.mat'), 'dvt_s', 'std_s', 'dvt_w', 'std_w', 'dvt_sl', 'std_sl');
        fprintf('Saved trial structs to %s\n', fullfile(outdir, 'R060721_trial_structs.mat'));
    catch ME
        warning(ME.identifier, '%s', ME.message);
    end

    %% Build the same numeric inputs used by the pipeline and compute MI
    fprintf('\nBuilding numeric inputs for single-run sleep/wake and comparison...\n');
    try
        [sleep_data, sleep_stim, wake_data, wake_stim, comp_data_A, comp_data_B, comp_stim] = build_numeric_inputs(dvt_s, std_s, dvt_w, std_w, dvt_sl, std_sl, baseline, chan);
    catch ME
        error('Failed to build numeric inputs: %s', ME.message);
    end

    fprintf('Computing MI for single-run sleep (cnm_MI_stimtime on numeric matrix)...\n');
    try
        [MI_sleep, sig_sleep] = cnm_MI_stimtime(sleep_data, sleep_stim, 0);
    catch ME
        error('Unable to run cnm_MI_stimtime for sleep: %s', ME.message);
    end

    fprintf('Computing MI for single-run wake (cnm_MI_stimtime on numeric matrix)...\n');
    try
        [MI_wake, sig_wake] = cnm_MI_stimtime(wake_data, wake_stim, 0);
    catch ME
        error('Unable to run cnm_MI_stimtime for wake: %s', ME.message);
    end

    fprintf('Computing matched-trial single-run MIs using the same trim as comparison...\n');
    try
        [MI_sleep_matched, sig_sleep_matched] = cnm_MI_stimtime(comp_data_B, comp_stim, 0);
        [MI_wake_matched, sig_wake_matched] = cnm_MI_stimtime(comp_data_A, comp_stim, 0);
    catch ME
        error('Unable to run matched-trial cnm_MI_stimtime: %s', ME.message);
    end

    save(fullfile(outdir,'R060721_single_MI_E14.mat'),'MI_sleep','MI_wake');
    fprintf('Saved single-run MI to %s\n', fullfile(outdir,'R060721_single_MI_E14.mat'));

    %% Compute comparison MI using comparison routine
    fprintf('\nComputing comparison MI (cnm_MI_stimtime_perm_diff on numeric matrices)...\n');
    try
        [MI_diff, MI_1, MI_2, sigMask_cmp] = cnm_MI_stimtime_perm_diff(comp_data_A, comp_data_B, comp_stim, 0, 0.05);
    catch ME
        error('Unable to run cnm_MI_stimtime_perm_diff: %s', ME.message);
    end
    MI_cmp = struct('MI_diff', MI_diff, 'MI_1', MI_1, 'MI_2', MI_2, 'sigMask', sigMask_cmp);
    save(fullfile(outdir,'R060721_comparison_MI_E14.mat'),'MI_cmp');
    fprintf('Saved comparison MI to %s\n', fullfile(outdir,'R060721_comparison_MI_E14.mat'));

    %% Helper to extract numeric vector for electrode E14
    % Note: local function implementation moved to file end to comply with
    % MATLAB script/function ordering rules.

    % Extract numeric vectors
    v_sleep = extract_E1(MI_sleep);
    v_wake  = extract_E1(MI_wake);

    % Comparison diff: try MI_cmp.MI_diff or MI_cmp.MI_2 - MI_cmp.MI_1
    v_cmp = [];
    if isstruct(MI_cmp)
        if isfield(MI_cmp,'MI_diff'), v_cmp = extract_E1(MI_cmp.MI_diff); end
        if isempty(v_cmp) && isfield(MI_cmp,'MI_2') && isfield(MI_cmp,'MI_1')
            a = MI_cmp.MI_2; b = MI_cmp.MI_1;
            if isstruct(a) && isfield(a,'E1') && isstruct(b) && isfield(b,'E1')
                v_cmp = a.E1 - b.E1;
            end
        end
        if isempty(v_cmp)
            % try direct extraction
            v_cmp = extract_E1(MI_cmp);
        end
    else
        v_cmp = extract_E1(MI_cmp);
    end

    %% Legacy CoI comparisons (single sleep file and comparison file)
    legacy_sleep_path = fullfile(repoRoot,'Data','Drosophila_CoI','R060721_sleep_BSLEEP','R060721_E14_permuted_E14.mat');
    legacy_comp_path  = fullfile(repoRoot,'Data','Drosophila_CoI','Comparisons','R060721_wake_sleep_BSLEEP','R060721_E14_permuted_E14.mat');

    legacy_sleep_MI1 = [];
    legacy_sleep_MI2 = [];
    legacy_comp_MI1 = [];
    legacy_comp_MI2 = [];

    if exist(legacy_sleep_path,'file')
        try
            Ls = load(legacy_sleep_path);
            leaf = Ls.CoI.R060721_E14_permuted_E14.CoI.R060721.BSLEEP;
            legacy_sleep_MI1 = leaf.MI1.E14_E14;
            legacy_sleep_MI2 = leaf.MI2.E14_E14;
        catch ME
            warning(ME.identifier, '%s', ME.message);
        end
    else
        warning('Legacy sleep CoI file not found at %s', legacy_sleep_path);
    end

    if exist(legacy_comp_path,'file')
        try
            Lc = load(legacy_comp_path);
            leaf = Lc.CoI.R060721_E14_permuted_E14.CoI.R060721.BSLEEP;
            legacy_comp_MI1 = leaf.MI1_diff.E14_E14;
            legacy_comp_MI2 = leaf.MI2_diff.E14_E14;
        catch ME
            warning(ME.identifier, '%s', ME.message);
        end
    else
        warning('Legacy comparison CoI file not found at %s', legacy_comp_path);
    end

    %% Compute stats and comparisons
    fprintf('\n--- Summary stats ---\n');

    print_stats('single_sleep_MI', v_sleep);
    print_stats('single_wake_MI', v_wake);
    print_stats('matched_sleep_MI', extract_E1(MI_sleep_matched));
    print_stats('matched_wake_MI', extract_E1(MI_wake_matched));
    if ~isempty(v_sleep) && ~isempty(v_wake) && isequal(size(v_sleep), size(v_wake))
        sd = double(v_sleep) - double(v_wake);
        fprintf('single-derived diff: mean=%g max=%g RMSE=%g\n', mean(sd(:),'omitnan'), max(abs(sd(:))), sqrt(mean(sd(:).^2,'omitnan')));
    else
        fprintf('single-derived diff: unavailable (missing or shape mismatch)\n');
    end

    matched_sleep = extract_E1(MI_sleep_matched);
    matched_wake = extract_E1(MI_wake_matched);
    if ~isempty(matched_sleep) && ~isempty(matched_wake) && isequal(size(matched_sleep), size(matched_wake))
        sdm = double(matched_sleep) - double(matched_wake);
        fprintf('matched-trial diff: mean=%g max=%g RMSE=%g\n', mean(sdm(:),'omitnan'), max(abs(sdm(:))), sqrt(mean(sdm(:).^2,'omitnan')));
    else
        fprintf('matched-trial diff: unavailable (missing or shape mismatch)\n');
    end

    % Compute co-information based MI1_diff on the matched trials with kperm=20.
    fprintf('\nComputing co-information based MI1_diff on matched trials (trimming BEFORE baseline, with kperm=20)...\n');
    try
        % First compute deterministic MI and coI for activity A and B
        [MI1_A, MI2_A, coI_A, jointMI_A, ~] = cnm_coI_stimtimetime_within_channel(comp_data_A, comp_data_A, comp_stim, 0);
        [MI1_B, MI2_B, coI_B, jointMI_B, ~] = cnm_coI_stimtimetime_within_channel(comp_data_B, comp_data_B, comp_stim, 0);

        MI1_diff_coi = MI1_A - MI1_B;
        MI2_diff_coi = MI2_A - MI2_B;
        diffCoI_matched = coI_B - coI_A;

        % Now run permutations locally with progress printing
        kperm_local = 20;
        nSamples = size(comp_data_A, 2);
        diffCoI_perm = nan(kperm_local, nSamples, nSamples);
        fprintf('Running %d permutations: [', kperm_local);
        for permIdx = 1:kperm_local
            fprintf('#');
            randorder = randperm(size(comp_data_A,1));
            [~, ~, coI_perm_A, ~, ~] = cnm_coI_stimtimetime_within_channel(comp_data_A, comp_data_A, comp_stim(randorder), 0);
            [~, ~, coI_perm_B, ~, ~] = cnm_coI_stimtimetime_within_channel(comp_data_B, comp_data_B, comp_stim(randorder), 0);
            diffCoI_perm(permIdx, :, :) = coI_perm_B - coI_perm_A;
        end
        fprintf(']\n');

        sigMask_coi = [];
    catch ME
        warning(ME.identifier, 'co-information matched-trial computation failed: %s', ME.message);
        MI1_diff_coi = [];
        sigMask_coi = [];
    end

    v_coi_MI1 = [];
    if ~isempty(MI1_diff_coi)
        v_coi_MI1 = extract_E1(MI1_diff_coi);
    end
    print_stats('coI_matched_MI1_diff', v_coi_MI1);
    if ~isempty(v_coi_MI1) && ~isempty(legacy_comp_MI1) && isequal(size(v_coi_MI1), size(legacy_comp_MI1))
        dco = double(v_coi_MI1) - double(legacy_comp_MI1);
        fprintf('coI matched MI1_diff - legacy MI1_diff: mean=%g max=%g RMSE=%g\n', mean(dco(:),'omitnan'), max(abs(dco(:))), sqrt(mean(dco(:).^2,'omitnan')));
    else
        fprintf('coI matched vs legacy MI1_diff comparison: unavailable (missing or shape mismatch)\n');
    end

    % Now reproduce legacy ordering: trim trials BEFORE baseline-normalisation
    fprintf('\nRecomputing co-information MI1_diff with trimming BEFORE baseline-normalisation (legacy order)...\n');
    try
        % determine trialnum like the legacy driver does
        trialnum_pre = min(std_w.trials, std_sl.trials);
        % trim raw structs
        dvt_w_pre = dvt_w; std_w_pre = std_w; dvt_sl_pre = dvt_sl; std_sl_pre = std_sl;
        dvt_w_pre.data = dvt_w_pre.data(:,:,1:trialnum_pre); dvt_w_pre.trials = trialnum_pre;
        std_w_pre.data = std_w_pre.data(:,:,1:trialnum_pre); std_w_pre.trials = trialnum_pre;
        dvt_sl_pre.data = dvt_sl_pre.data(:,:,1:trialnum_pre); dvt_sl_pre.trials = trialnum_pre;
        std_sl_pre.data = std_sl_pre.data(:,:,1:trialnum_pre); std_sl_pre.trials = trialnum_pre;

        % baseline-normalise the trimmed structs via the same prep_condition
        [comp_data_A_pre, comp_stim_pre] = prep_condition(dvt_w_pre, std_w_pre, baseline, chan);
        [comp_data_B_pre, ~] = prep_condition(dvt_sl_pre, std_sl_pre, baseline, chan);

        % run the within-channel coI perm-diff with kperm=20 and a progress bar
        nperms = 20;
        fprintf('Using kperm=%d\n', nperms);
        [MI1_A_pre, MI2_A_pre, coI_A_pre, jointMI_A_pre, ~] = cnm_coI_stimtimetime_within_channel(comp_data_A_pre, comp_data_A_pre, comp_stim_pre, 0);
        [MI1_B_pre, MI2_B_pre, coI_B_pre, jointMI_B_pre, ~] = cnm_coI_stimtimetime_within_channel(comp_data_B_pre, comp_data_B_pre, comp_stim_pre, 0);
        MI1_diff_coi_pre = MI1_A_pre - MI1_B_pre;
        MI2_diff_coi_pre = MI2_A_pre - MI2_B_pre;
        diffCoI_pre = coI_B_pre - coI_A_pre;
        diffJointMI_pre = jointMI_A_pre - jointMI_B_pre;
        nSamples_pre = size(comp_data_A_pre, 2);
        diffCoI_perm_pre = nan(nperms, nSamples_pre, nSamples_pre);
        fprintf('Running %d permutations: [', nperms);
        for permIdx = 1:nperms
            fprintf('#');
            randorder = randperm(size(comp_data_A_pre,1));
            [~, ~, coI_perm_A_pre, ~, ~] = cnm_coI_stimtimetime_within_channel(comp_data_A_pre, comp_data_A_pre, comp_stim_pre(randorder), 0);
            [~, ~, coI_perm_B_pre, ~, ~] = cnm_coI_stimtimetime_within_channel(comp_data_B_pre, comp_data_B_pre, comp_stim_pre(randorder), 0);
            diffCoI_perm_pre(permIdx, :, :) = coI_perm_B_pre - coI_perm_A_pre;
        end
        fprintf(']\n');
        sigMask_pre = [];
    catch ME
        warning(ME.identifier, 'Pretrim coI computation failed: %s', ME.message);
        MI1_diff_coi_pre = [];
    end
    v_coi_MI1_pre = [];
    if ~isempty(MI1_diff_coi_pre)
        v_coi_MI1_pre = extract_E1(MI1_diff_coi_pre);
    end
    print_stats('coI_pretrim_MI1_diff', v_coi_MI1_pre);
    if ~isempty(v_coi_MI1_pre) && ~isempty(legacy_comp_MI1) && isequal(size(v_coi_MI1_pre), size(legacy_comp_MI1))
        dpre = double(v_coi_MI1_pre) - double(legacy_comp_MI1);
        fprintf('coI pretrim MI1_diff - legacy MI1_diff: mean=%g max=%g RMSE=%g\n', mean(dpre(:),'omitnan'), max(abs(dpre(:))), sqrt(mean(dpre(:).^2,'omitnan')));
    else
        fprintf('coI pretrim vs legacy MI1_diff comparison: unavailable (missing or shape mismatch)\n');
    end

    print_stats('comparison_diff', v_cmp);
    print_stats('legacy_sleep_MI1', legacy_sleep_MI1);
    print_stats('legacy_sleep_MI2', legacy_sleep_MI2);
    print_stats('legacy_comp_MI1_diff', legacy_comp_MI1);
    print_stats('legacy_comp_MI2_diff', legacy_comp_MI2);

    if ~isempty(v_sleep) && ~isempty(legacy_sleep_MI1) && isequal(size(v_sleep), size(legacy_sleep_MI1))
        ds = double(v_sleep) - double(legacy_sleep_MI1);
        fprintf('single sleep - legacy sleep MI1: mean=%g max=%g RMSE=%g\n', mean(ds(:),'omitnan'), max(abs(ds(:))), sqrt(mean(ds(:).^2,'omitnan')));
    else
        fprintf('single sleep vs legacy sleep MI1 comparison: unavailable (missing or shape mismatch)\n');
    end

    if ~isempty(v_wake) && ~isempty(legacy_sleep_MI2) && isequal(size(v_wake), size(legacy_sleep_MI2))
        dw = double(v_wake) - double(legacy_sleep_MI2);
        fprintf('single wake - legacy sleep MI2: mean=%g max=%g RMSE=%g\n', mean(dw(:),'omitnan'), max(abs(dw(:))), sqrt(mean(dw(:).^2,'omitnan')));
    else
        fprintf('single wake vs legacy sleep MI2 comparison: unavailable (missing or shape mismatch)\n');
    end

    if ~isempty(legacy_comp_MI1)
        if ~isempty(v_cmp) && isequal(size(v_cmp), size(legacy_comp_MI1))
            dcl = double(v_cmp) - double(legacy_comp_MI1);
            fprintf('comp diff - legacy MI1_diff: mean=%g max=%g RMSE=%g\n', mean(dcl(:),'omitnan'), max(abs(dcl(:))), sqrt(mean(dcl(:).^2,'omitnan')));
        else
            fprintf('comp vs legacy MI1_diff comparison: unavailable (missing or shape mismatch)\n');
        end

        if ~isempty(matched_sleep) && ~isempty(matched_wake) && isequal(size(matched_sleep), size(matched_wake)) && isequal(size(matched_sleep), size(legacy_comp_MI1))
            matched_diff = double(matched_sleep) - double(matched_wake);
            dman = matched_diff - double(legacy_comp_MI1);
            fprintf('matched diff - legacy MI1_diff: mean=%g max=%g RMSE=%g\n', mean(dman(:),'omitnan'), max(abs(dman(:))), sqrt(mean(dman(:).^2,'omitnan')));
        else
            fprintf('matched diff vs legacy MI1_diff comparison: unavailable (missing or shape mismatch)\n');
        end
    else
        fprintf('legacy comparison MI1_diff unavailable\n');
    end

    if ~isempty(v_cmp) && ~isempty(v_sleep) && ~isempty(v_wake)
        try
            sd = double(v_sleep) - double(v_wake);
            if isequal(size(v_cmp), size(sd))
                dcs = double(v_cmp) - double(sd);
                fprintf('comp - single-derived: mean=%g max=%g RMSE=%g\n', mean(dcs(:),'omitnan'), max(abs(dcs(:))), sqrt(mean(dcs(:).^2,'omitnan')));
            else
                fprintf('comp vs single-derived shapes differ: comp=%s single=%s\n', mat2str(size(v_cmp)), mat2str(size(sd)));
            end
        catch
        end
    end

    % Save everything
    try
        save(fullfile(outdir,'R060721_MI_all.mat'),'v_sleep','v_wake','matched_sleep','matched_wake','v_cmp','legacy_sleep_MI1','legacy_sleep_MI2','legacy_comp_MI1','legacy_comp_MI2');
        fprintf('Saved combined results to %s\n', fullfile(outdir,'R060721_MI_all.mat'));
    catch
        warning('Failed saving combined results');
    end

    fprintf('Done.\n');
catch ME
    fprintf('Script failed: %s\n', ME.message);
    for k=1:length(ME.stack)
        s = ME.stack(k);
        fprintf('  %s:%d\n', s.file, s.line);
    end
    rethrow(ME);
end

end % end of main function

function v = extract_E1(x)
% Recursively find a numeric vector for channel E1 inside various possible
% structures returned by the pipeline functions.
    v = [];
    if isempty(x), return; end
    if isstruct(x)
        if isfield(x,'E1')
            v = x.E1; return;
        end
        if isfield(x,'MI_2') && isfield(x,'MI_1')
            a = x.MI_2; b = x.MI_1;
            if isstruct(a) && isfield(a,'E1') && isstruct(b) && isfield(b,'E1')
                v = a.E1 - b.E1; return;
            end
        end
        fn = fieldnames(x);
        for ii = 1:length(fn)
            try
                vv = extract_E1(x.(fn{ii}));
                if ~isempty(vv), v = vv; return; end
            catch
            end
        end
    elseif isnumeric(x)
        v = x; return;
    else
        try
            fn = fieldnames(x);
            for ii = 1:length(fn)
                vv = extract_E1(x.(fn{ii}));
                if ~isempty(vv), v = vv; return; end
            end
        catch
        end
    end
end

function print_stats(name, arr)
    if isempty(arr)
        fprintf('%s: not found\n', name);
        return;
    end
    arrnum = double(arr(:));
    fprintf('%s: size=%s mean=%g max=%g\n', name, mat2str(size(arr)), mean(arrnum,'omitnan'), max(arrnum));
end

function [sleep_data, sleep_stim, wake_data, wake_stim, comp_data_A, comp_data_B, comp_stim] = build_numeric_inputs(dvt_s, std_s, dvt_w, std_w, dvt_sl, std_sl, baseline, chan)
% Reproduce the main pipeline's data preparation so cnm_MI_stimtime sees
% numeric trial x time matrices and a simple 2-class stimulus vector.

    [sleep_data, sleep_stim] = prep_condition(dvt_s, std_s, baseline, chan);
    [wake_data, wake_stim] = prep_condition(dvt_w, std_w, baseline, chan);
    [comp_data_A, comp_stim] = prep_condition(dvt_w, std_w, baseline, chan);
    [comp_data_B, ~] = prep_condition(dvt_sl, std_sl, baseline, chan);

    % Match the comparison pipeline: trim both conditions to the same
    % number of trials before computing MI differences.
    comp_trialnum = min(size(comp_data_A, 1), size(comp_data_B, 1));
    comp_data_A = comp_data_A(1:comp_trialnum, :);
    comp_data_B = comp_data_B(1:comp_trialnum, :);
    comp_stim = comp_stim(1:comp_trialnum);
    fprintf('Comparison trialnum after trim: %d\n', comp_trialnum);
end

function [data_mat, stim] = prep_condition(dvt, std, baseline, chan)
    if ~isfield(dvt,'data') || ~isfield(std,'data')
        error('Expected dvt/std structs with data fields');
    end

    % Match the preprocessing in main_MI_ERP.m / main_MI_ERP_comparison.m
    data_dvt = permute(dvt.data,[1 3 2]);
    mean_data_dvt = squeeze(mean(data_dvt(:,:,baseline),3));
    dvt.data = data_dvt - repmat(mean_data_dvt,[1 1 size(dvt.data,2)]);
    bb1_dev_bl = permute(dvt.data,[1 3 2]);

    data_std = permute(std.data,[1 3 2]);
    mean_data_std = squeeze(mean(data_std(:,:,baseline),3));
    std.data = data_std - repmat(mean_data_std,[1 1 size(std.data,2)]);
    bb1_std_bl = permute(std.data,[1 3 2]);

    bb_dev = permute(bb1_dev_bl,[3 1 2]);
    bb_std = permute(bb1_std_bl,[3 1 2]);

    % Use the selected channel.
    data_mat = [squeeze(bb_dev(:,chan,:)); squeeze(bb_std(:,chan,:))];
    stim = [zeros(size(bb_dev,1),1); ones(size(bb_std,1),1)];
end
