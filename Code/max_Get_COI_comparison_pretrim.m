%% GET CO-INFORMATION
% Pretrim version of max_Get_COI_comparison.
% Trials are randomly matched before baseline normalisation so the more
% numerous activity is downsampled without taking the first N trials.

function [MI_stat] = max_Get_COI_comparison_pretrim(electrodes, basefold, datatype,...
    subject, condition, deviant_group_number_1, ...
    standard_group_number_1, deviant_group_number_2, standard_group_number_2,...
    corrected, srate, baseline, kperm, threshold)

rng(0, 'twister');

overVar_file = strcat(basefold, datatype, '/', subject, '_', condition, '.mat');
[dvt_1, std_1] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_1, standard_group_number_1,...
    corrected, srate);
[dvt_2, std_2] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_2, standard_group_number_2,...
    corrected, srate);

patname = char(subject);

%% MATCH TRIALS BEFORE BASELINE NORMALISATION
[dvt_1, std_1, dvt_2, std_2, trialnum] = pretrim_match_trials(...
    dvt_1, std_1, dvt_2, std_2, ...
    length(deviant_group_number_1), length(standard_group_number_1), ...
    length(deviant_group_number_2), length(standard_group_number_2));

fprintf('Pretrim CoI for %s %s: using %d matched trials\n', patname, condition, trialnum);

%% BASELINE NORMALISE AFTER MATCHING
data_dvt_1 = permute(dvt_1.data,[1 3 2]); mean_data_dvt_1 = squeeze(mean(data_dvt_1(:,:,baseline),3));
dvt_1.data = data_dvt_1 - repmat(mean_data_dvt_1,[1 1 size(dvt_1.data,2)]);
bb1_dev_bl_1 = permute(dvt_1.data,[1 3 2]);

data_std_1 = permute(std_1.data,[1 3 2]); mean_data_std_1 = squeeze(mean(data_std_1(:,:,baseline),3));
std_1.data = data_std_1 - repmat(mean_data_std_1,[1 1 size(std_1.data,2)]);
bb1_std_bl_1 = permute(std_1.data,[1 3 2]);

dvt_1_data_minus_mean = permute(bb1_dev_bl_1,[3 1 2]);
std_1_data_minus_mean = permute(bb1_std_bl_1,[3 1 2]);

data_dvt_2 = permute(dvt_2.data,[1 3 2]); mean_data_dvt_2 = squeeze(mean(data_dvt_2(:,:,baseline),3));
dvt_2.data = data_dvt_2 - repmat(mean_data_dvt_2,[1 1 size(dvt_2.data,2)]);
bb1_dev_bl_2 = permute(dvt_2.data,[1 3 2]);

data_std_2 = permute(std_2.data,[1 3 2]); mean_data_std_2 = squeeze(mean(data_std_2(:,:,baseline),3));
std_2.data = data_std_2 - repmat(mean_data_std_2,[1 1 size(std_2.data,2)]);
bb1_std_bl_2 = permute(std_2.data,[1 3 2]);

dvt_2_data_minus_mean = permute(bb1_dev_bl_2,[3 1 2]);
std_2_data_minus_mean = permute(bb1_std_bl_2,[3 1 2]);

%% GET THE RIGHT CHANNELS
if isempty(std_1.chanlocs) == 1
    all_elecs = {};
    for e_index = 1:std_1.nbchan
        name = strcat('E', num2str(e_index));
        all_elecs(e_index) = {name};
    end
    All_elecs = all_elecs';
else
    All_elecs = std_1.chanlocs.labels;
end

channel1 = extractBefore(electrodes,'_permuted_');
channel2 = extractAfter(electrodes,'_permuted_');
chnum1 = extractBefore(electrodes,'_permuted_');
chnum2 = extractAfter(electrodes,'_permuted_');
channel1 = strcmp(channel1, All_elecs');
channel2 = strcmp(channel2, All_elecs');
comb_chan = char(strcat(chnum1,'_',chnum2));
comb_chan = replace(comb_chan,'-','_');

dvt_1_data_ERP_sel_chan1 = squeeze(dvt_1_data_minus_mean(:, channel1,:));
std_1_data_ERP_sel_chan1 = squeeze(std_1_data_minus_mean(:,channel1,:));
dvt_1_data_ERP_sel_chan2 = squeeze(dvt_1_data_minus_mean(:, channel2, :));
std_1_data_ERP_sel_chan2 = squeeze(std_1_data_minus_mean(:,channel2,:));

dvt_2_data_ERP_sel_chan1 = squeeze(dvt_2_data_minus_mean(:, channel1,:));
std_2_data_ERP_sel_chan1 = squeeze(std_2_data_minus_mean(:,channel1,:));
dvt_2_data_ERP_sel_chan2 = squeeze(dvt_2_data_minus_mean(:, channel2, :));
std_2_data_ERP_sel_chan2 = squeeze(std_2_data_minus_mean(:,channel2,:));

%% COMPUTE CO-INFORMATION WITH GCMI
if strcmp(chnum1,chnum2) == 1
    [MI1_diff, MI2_diff, diffCoI, diffJointMI, sigMask] = cnm_coI_stimtimetime_perm_diff_within_channel(...
     [dvt_1_data_ERP_sel_chan1; std_1_data_ERP_sel_chan1],...
     [dvt_1_data_ERP_sel_chan2; std_1_data_ERP_sel_chan2],...
     [dvt_2_data_ERP_sel_chan1; std_2_data_ERP_sel_chan1],...
     [dvt_2_data_ERP_sel_chan2; std_2_data_ERP_sel_chan2],...
     [zeros(1, size(dvt_1_data_minus_mean,1)), ones(1, size(dvt_1_data_minus_mean,1))]', kperm, threshold);
else
    [MI1_diff, MI2_diff, diffCoI, diffJointMI, sigMask] = cnm_coI_stimtimetime_perm_diff(...
     [dvt_1_data_ERP_sel_chan1; std_1_data_ERP_sel_chan1],...
     [dvt_1_data_ERP_sel_chan2; std_1_data_ERP_sel_chan2],...
     [dvt_2_data_ERP_sel_chan1; std_2_data_ERP_sel_chan1],...
     [dvt_2_data_ERP_sel_chan2; std_2_data_ERP_sel_chan2],...
     [zeros(1, size(dvt_1_data_minus_mean,1)), ones(1, size(dvt_1_data_minus_mean,1))]', kperm, threshold);
end

MI_stat.CoI.(patname).(char(condition)).diffCoI.(comb_chan) = diffCoI;
MI_stat.CoI.(patname).electrode.(char(condition)).(comb_chan) = electrodes;
MI_stat.CoI.(patname).(char(condition)).sigMask.(comb_chan) = sigMask;
MI_stat.CoI.(patname).(char(condition)).MI1_diff.(comb_chan) = MI1_diff;
MI_stat.CoI.(patname).(char(condition)).MI2_diff.(comb_chan) = MI2_diff;
MI_stat.CoI.(patname).(char(condition)).jointJointMI.(comb_chan) = diffJointMI;
end

function [dvt_1, std_1, dvt_2, std_2, trialnum] = pretrim_match_trials(dvt_1, std_1, dvt_2, std_2, n_dev_1, n_std_1, n_dev_2, n_std_2)
% Randomly match trials before baseline normalisation.
% If the deviant/standard groups are color-balanced pairs, preserve that
% balance by sampling the same number from each sub-group.
if n_dev_1 == 2 && n_dev_2 == 2 && n_std_1 == 2 && n_std_2 == 2
    per_group_target = floor(min([size(dvt_1.data,3)/n_dev_1, size(std_1.data,3)/n_std_1, ...
        size(dvt_2.data,3)/n_dev_2, size(std_2.data,3)/n_std_2]));

    if per_group_target < 1
        error('Not enough trials to perform subtype-balanced matching.');
    end

    dvt_1.data = sample_equal_from_concatenated_groups(dvt_1.data, n_dev_1, per_group_target);
    std_1.data = sample_equal_from_concatenated_groups(std_1.data, n_std_1, per_group_target);
    dvt_2.data = sample_equal_from_concatenated_groups(dvt_2.data, n_dev_2, per_group_target);
    std_2.data = sample_equal_from_concatenated_groups(std_2.data, n_std_2, per_group_target);

    dvt_1.trials = size(dvt_1.data,3);
    std_1.trials = size(std_1.data,3);
    dvt_2.trials = size(dvt_2.data,3);
    std_2.trials = size(std_2.data,3);
    trialnum = std_1.trials;
else
    [dvt_1, std_1, dvt_2, std_2, trialnum] = random_match_trials_total(dvt_1, std_1, dvt_2, std_2);
end
end

function data_out = sample_equal_from_concatenated_groups(data_in, n_groups, per_group_target)
% Assumes trials are concatenated by group in contiguous blocks.
total_trials = size(data_in,3);
if mod(total_trials, n_groups) ~= 0
    error('Trial count (%d) is not divisible by n_groups (%d).', total_trials, n_groups);
end

trials_per_group = total_trials / n_groups;
all_idx = [];
for g = 1:n_groups
    first_idx = (g - 1) * trials_per_group + 1;
    last_idx = g * trials_per_group;
    group_idx = first_idx:last_idx;
    pick_local = randperm(trials_per_group, per_group_target);
    all_idx = [all_idx, group_idx(pick_local)]; %#ok<AGROW>
end

all_idx = all_idx(randperm(length(all_idx)));
data_out = data_in(:,:,all_idx);
end

function [dvt_1, std_1, dvt_2, std_2, trialnum] = random_match_trials_total(dvt_1, std_1, dvt_2, std_2)
% Fallback: random total-trial matching without subtype balancing.
trialnum = min([std_1.trials, std_2.trials, dvt_1.trials, dvt_2.trials]);

idx_std_1 = randperm(std_1.trials, trialnum);
idx_dvt_1 = randperm(dvt_1.trials, trialnum);
idx_std_2 = randperm(std_2.trials, trialnum);
idx_dvt_2 = randperm(dvt_2.trials, trialnum);

std_1.data = std_1.data(:,:,idx_std_1);
dvt_1.data = dvt_1.data(:,:,idx_dvt_1);
std_2.data = std_2.data(:,:,idx_std_2);
dvt_2.data = dvt_2.data(:,:,idx_dvt_2);

std_1.trials = trialnum;
dvt_1.trials = trialnum;
std_2.trials = trialnum;
dvt_2.trials = trialnum;
end