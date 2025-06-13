function correlation_main(task_id)
%%CORRELATION_MAIN Entry point for electrode pair correlation analysis
% This function is called by array jobs with task_id from 1-840
% Each task_id corresponds to a specific fly and electrode pair

tic
% Define parameters
USING_HPC = 2; % Set to 2 for Apocrita

%% Step 1: Add paths and get parameters
if USING_HPC == 2 % Add paths for QMUL Apocrita
    addpath(genpath('/data/SBBS-PIDProject/Maxime/CNM'));
    addpath(genpath('/data/SBBS-PIDProject/Maxime/GCMI_master'));
    addpath(genpath('/data/SBBS-PIDProject/Maxime/Prediction_Errors'));
end

% Get parameters for the analysis
[basefold, datatype, all_con, condition, subject, participants, EoI,...
    srate, activity_tag, deviant_group_number, standard_group_number, corrected, stim_onset, baseline,...
    start_cut_off, end_cut_off, kperm] = Max_get_param(USING_HPC, 0);

%% Step 2: Map task_id to fly and electrode pair
% With 8 flies, 15 electrodes per fly (105 pairs per fly)
% Total tasks = 8 * 105 = 840
num_electrodes = 15;
pairs_per_fly = nchoosek(num_electrodes, 2); % 105 pairs per fly

% Map task_id to fly and electrode pair
[fly_idx, electrode_pair] = map_task_id_to_fly_and_pair(task_id, pairs_per_fly);

% Get fly name from participants list
fly_name = participants{fly_idx};
e1 = electrode_pair(1);
e2 = electrode_pair(2);

fprintf('Processing fly %s (%d of 8), electrodes E%d-E%d (task_id: %d)\n', fly_name, fly_idx, e1, e2, task_id);

%% Step 3: Run correlation analysis for this electrode pair
% Load data
data_file = fullfile(basefold, datatype, [fly_name, '_', condition, '.mat']);
if ~exist(data_file, 'file')
    error('Data file not found: %s', data_file);
end

% Load deviant and standard trials
[dvt, std] = load_trials_from_group_hyper(data_file, deviant_group_number, standard_group_number, corrected, srate);

% Extract data for selected electrodes
e1_dev = squeeze(dvt.data(e1, :, :))';  % Transpose to get Trials x Time
e2_dev = squeeze(dvt.data(e2, :, :))';
e1_std = squeeze(std.data(e1, :, :))';
e2_std = squeeze(std.data(e2, :, :))';

% Set parameters for correlation analysis
max_lag = 20; % Maximum lag in samples

% Calculate cross-correlation
[corr_values, lags, max_corr, optimal_lag, trial_corrs] = ...
    electrode_diff_correlation(e1_dev, e1_std, e2_dev, e2_std, max_lag, srate, false);

%% Step 4: Save results
% Create results structure
results = struct();
results.electrode_pair = [e1, e2];
results.corr_values = corr_values;
results.lags = lags;
results.max_corr = max_corr;
results.optimal_lag = optimal_lag;
results.trial_corrs = trial_corrs;
results.fly_name = fly_name;
results.condition = condition;
results.srate = srate;

% Set output directory
if USING_HPC == 2
    results_dir = sprintf('/data/SBBS-PIDProject/Maxime/Drosophila_Results/Correlation_Analysis/%s_%s', ...
        fly_name, activity_tag);
else
    results_dir = fullfile(basefold, 'Correlation_Analysis', sprintf('%s_%s', fly_name, activity_tag));
end

% Create output directory if it doesn't exist
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Save results
filename = sprintf('%s_%s_E%d_E%d_correlation.mat', fly_name, condition, e1, e2);
save_path = fullfile(results_dir, filename);
save(save_path, 'results');

fprintf('Completed analysis for fly %s, E%d-E%d. Results saved to %s\n', fly_name, e1, e2, save_path);
toc
end

function [fly_idx, pair] = map_task_id_to_fly_and_pair(task_id, pairs_per_fly)
% MAP_TASK_ID_TO_FLY_AND_PAIR Maps task_id (1-840) to a fly and electrode pair
% Input:
%   task_id      - Task ID (1-840)
%   pairs_per_fly - Number of electrode pairs per fly (105)
% Output:
%   fly_idx      - Index of the fly (1-8)
%   pair         - Electrode pair [e1, e2]

if task_id < 1 || task_id > 8 * pairs_per_fly
    error('Invalid task_id: %d. Must be between 1 and %d', task_id, 8 * pairs_per_fly);
end

% Calculate fly index (1-8)
fly_idx = ceil(task_id / pairs_per_fly);

% Calculate pair index within this fly (1-105)
pair_idx = task_id - (fly_idx - 1) * pairs_per_fly;
if pair_idx == 0
    pair_idx = pairs_per_fly;
    fly_idx = fly_idx - 1;
end

% Map pair_idx to electrode pair
num_electrodes = 15;
pair_count = 0;
pair = [0, 0];

for e1 = 1:num_electrodes
    for e2 = (e1+1):num_electrodes
        pair_count = pair_count + 1;
        if pair_count == pair_idx
            pair = [e1, e2];
            return;
        end
    end
end
end