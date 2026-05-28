%% ESTIMATES MI AND GETS ELECTRODES OF INTEREST
% Copy of main_MI_ERP_comparison with the preprocessing order made explicit:
% trim trials first, then baseline-normalise, then estimate MI.

function main_MI_ERP_comparison_pretrim(task_id)
USING_HPC = 2; % 0 for local, 1 for Cambridge HPC, 2 for QMUL HPC
if USING_HPC == 1
    addpath(genpath('/home/mj649/rds/hpc-work/CNM'));
    addpath(genpath('/home/mj649/rds/hpc-work/GCMI_master'));
    addpath(genpath('/home/mj649/rds/hpc-work/Prediction_Errors'));
end

if USING_HPC == 2
    addpath(genpath('/data/SBBS-PIDProject/Maxime/CNM'));
    addpath(genpath('/data/SBBS-PIDProject/Maxime/GCMI_master'));
    addpath(genpath('/data/SBBS-PIDProject/Maxime/Prediction_Errors'));
end

THRESHOLD = 0.05;
rng(0, 'twister'); % Reproducible single random draw per run.
[basefold, datatype, all_con, condition, subject, participants, EoI, ...
    srate, activity_tag_1, deviant_group_number_1, standard_group_number_1, ...
    activity_tag_2, deviant_group_number_2, standard_group_number_2, corrected, ...
    stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_comparison_param(USING_HPC, 0);

for i = 1:length(participants)
    for con = 1:length(all_con)
        overVar_file = strcat(basefold, datatype, '/', participants(i), '_', all_con(con), '.mat');
        if exist(char(overVar_file)) == 0
            continue
        end

        [dvt_1, std_1] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_1, standard_group_number_1, corrected, srate);
        [dvt_2, std_2] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_2, standard_group_number_2, corrected, srate);

        [bb_dev_1, bb_std_1, bb_dev_2, bb_std_2, timing, trialnum] = prepare_trimmed_baselined_inputs(dvt_1, std_1, dvt_2, std_2, baseline);

        fprintf('Participant %s condition %s: using %d matched trials\n', char(participants(i)), char(all_con(con)), trialnum);

        elec_of_I = {};
        h = figure(con);
        for ch = 1:dvt_1.nbchan
            chnum = strcat('ch_', string(ch));
            dat_1.class1.BB = squeeze(bb_dev_1(:, ch, :));
            dat_1.class2.BB = squeeze(bb_std_1(:, ch, :));
            dat_2.class1.BB = squeeze(bb_dev_2(:, ch, :));
            dat_2.class2.BB = squeeze(bb_std_2(:, ch, :));

            [MI_diff, MI_1, MI_2, sigMask] = cnm_MI_stimtime_perm_diff( ...
                [dat_1.class1.BB; dat_1.class2.BB], ...
                [dat_2.class1.BB; dat_2.class2.BB], ...
                [zeros(1, size(bb_dev_1,1)), ones(1, size(bb_dev_1,1))]', ...
                kperm, THRESHOLD);

            chan_name = strcat('E', num2str(ch));
            MI_stat.MI_diff.(chan_name) = MI_diff;
            MI_stat.MI_1.(chan_name) = MI_1;
            MI_stat.MI_2.(chan_name) = MI_2;
            MI_stat.electrode = chan_name;
            MI_stat.sigMask.(chan_name) = sigMask;

            EoI_list = true;
            if EoI_list == true
                Electorodes.(chan_name) = chan_name;
            end
            hold off
        end

        MI_name = char(strcat(participants(i), '_pretrim_', (strcat(char(activity_tag_1), '_', char(activity_tag_2))), '_', char(all_con(con)), '_MI_data.mat'));
        save(fullfile(basefold, 'MI_Data', MI_name), 'MI_stat', '-mat');

        filename = fullfile(basefold, 'DataEoI', ['EoI_data_' datatype '_pretrim.mat']);
        if exist(filename, 'file')
            load(filename, 'EoI');
        else
            EoI = struct();
        end

        EoI.(char(participants(i))).(strcat(char(activity_tag_1), '_', char(activity_tag_2))).(char(all_con(con))) = Electorodes;
        EoI.(char(participants(i))).(strcat(char(activity_tag_1), '_', char(activity_tag_2))).(char(all_con(con))) = fieldnames(EoI.(char(participants(i))).(strcat(char(activity_tag_1), '_', char(activity_tag_2))).(char(all_con(con))));
        save(filename, 'EoI', '-mat');

        clear Electorodes
        clear MI_stat
        close all
    end
end
end

function [bb_dev_1, bb_std_1, bb_dev_2, bb_std_2, timing, trialnum] = prepare_trimmed_baselined_inputs(dvt_1, std_1, dvt_2, std_2, baseline)
% Match trials first using a single random subsample, then baseline-normalise.
% This avoids bias from always taking the first N trials.

if (std_1.trials < std_2.trials)
    trialnum = std_1.trials;
    idx_std_2 = randperm(std_2.trials, trialnum);
    idx_dvt_2 = randperm(dvt_2.trials, trialnum);
    std_2.data = std_2.data(:,:,idx_std_2);
    std_2.trials = trialnum;
    dvt_2.data = dvt_2.data(:,:,idx_dvt_2);
    dvt_2.trials = trialnum;
else
    trialnum = std_2.trials;
    idx_std_1 = randperm(std_1.trials, trialnum);
    idx_dvt_1 = randperm(dvt_1.trials, trialnum);
    std_1.data = std_1.data(:,:,idx_std_1);
    std_1.trials = trialnum;
    dvt_1.data = dvt_1.data(:,:,idx_dvt_1);
    dvt_1.trials = trialnum;
end

data_dvt_1 = permute(dvt_1.data,[1 3 2]);
mean_data_dvt_1 = squeeze(mean(data_dvt_1(:,:,baseline),3));
dvt_1.data = data_dvt_1 - repmat(mean_data_dvt_1,[1 1 size(dvt_1.data,2)]);
bb1_dev_bl_1 = permute(dvt_1.data,[1 3 2]);

data_std_1 = permute(std_1.data,[1 3 2]);
mean_data_std_1 = squeeze(mean(data_std_1(:,:,baseline),3));
std_1.data = data_std_1 - repmat(mean_data_std_1,[1 1 size(std_1.data,2)]);
bb1_std_bl_1 = permute(std_1.data,[1 3 2]);

bb_dev_1 = permute(bb1_dev_bl_1,[3 1 2]);
bb_std_1 = permute(bb1_std_bl_1,[3 1 2]);

data_dvt_2 = permute(dvt_2.data,[1 3 2]);
mean_data_dvt_2 = squeeze(mean(data_dvt_2(:,:,baseline),3));
dvt_2.data = data_dvt_2 - repmat(mean_data_dvt_2,[1 1 size(dvt_2.data,2)]);
bb1_dev_bl_2 = permute(dvt_2.data,[1 3 2]);

data_std_2 = permute(std_2.data,[1 3 2]);
mean_data_std_2 = squeeze(mean(data_std_2(:,:,baseline),3));
std_2.data = data_std_2 - repmat(mean_data_std_2,[1 1 size(std_2.data,2)]);
bb1_std_bl_2 = permute(std_2.data,[1 3 2]);

bb_dev_2 = permute(bb1_dev_bl_2,[3 1 2]);
bb_std_2 = permute(bb1_std_bl_2,[3 1 2]);

timing = dvt_1.times;
end