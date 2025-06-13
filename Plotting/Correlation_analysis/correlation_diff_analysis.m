% Get parameters from Max_get_param
[basefold, datatype, all_con, condition, subject, participants, EoI, ...
 srate, activity_tag, deviant_group_number, standard_group_number, corrected, ...
 stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_param(0, 0);

% Specify fly and electrode pair
fly_name = 'R060721';  
electrode_pair = [1, 3];  % Electrodes to analyze

% Load data
data_file = fullfile(basefold, datatype, [fly_name, '_', condition, '.mat']);
[dvt, std] = load_trials_from_group_hyper(data_file, deviant_group_number, standard_group_number, corrected, srate);

% Extract data for selected electrodes
e1_dev = squeeze(dvt.data(electrode_pair(1), :, :))';  % Transpose to get Trials x Time
e2_dev = squeeze(dvt.data(electrode_pair(2), :, :))';
e1_std = squeeze(std.data(electrode_pair(1), :, :))';
e2_std = squeeze(std.data(electrode_pair(2), :, :))';

% Calculate cross-correlation with confidence intervals
max_lag = 20; % Maximum lag in samples
[corr_values, lags, max_corr, optimal_lag, trial_corrs] = ...
    electrode_diff_correlation(e1_dev, e1_std, e2_dev, e2_std, max_lag, srate);

