%% CoI Comparison Analysis
% 
%% SETUP AND PARAMETERS
% Load CoI comparison data with significance masks from permutation testing
USING_HPC = 0;  % 0 = local machine, 1 = Cambridge HPC, 2 = QMUL HPC
get_elec = 0;   % 0 = don't get electrodes of interest (not needed for plotting)
times = -25:74; % Time vector in ms (baseline: -25 to 0, post-stim: 0 to 74)
stim_idx = 26;  % Index of stimulus onset (time = 0 ms)
n_flies = 8;    % Total number of flies in the dataset
cutoff = 8;     % Cutoff for peripheral/central division

%% LOAD PARAMETERS
% Get experimental parameters using the comparison parameter function
% This loads parameters for BOTH activity states (activity_tag_1 and activity_tag_2)
% Examples: wake vs sleep, beginning_sleep vs mid_sleep, etc.
[basefold, datatype, all_con, condition,subject, participants, EoI,...
srate, activity_tag_1, deviant_group_number_1, standard_group_number_1,...
activity_tag_2, deviant_group_number_2, standard_group_number_2, corrected,...
stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_comparison_param(USING_HPC, get_elec);

%% RUN ANALYSIS
[CoI_comparison_results] = compare_activity_conditions(participants, basefold, datatype, condition, cutoff, times, activity_tag_1, activity_tag_2);