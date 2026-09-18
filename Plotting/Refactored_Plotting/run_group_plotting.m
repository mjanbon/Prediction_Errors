% Script to run group-level CoI/MI/mask plotting

set_default_plotting();

% Parameters (adapt as needed)
USING_HPC = 0;
get_elec = 0;
cutoff = [5,11];
timing = -25:1:74; % 100 timepoints
times = 100;

plot_cfg = struct(...
    'xlimits', [-25 75], ...
    'ylimits_CoI', [-25 75], ...
    'ylimit_MI', [-0.005 0.1], ...
    'xticks_CoI', -25:25:75, ...
    'yticks_CoI', -25:25:75, ...
    'yticks_MI', [0:0.05:0.2], ...
    'x_labels', -25:25:75, ...
    'y_labels_CoI', -25:25:75, ...
    'y_labels_MI', [0:0.05:0.1], ...
    'climits', [-0.01 0.01], ...
    'climits_mask', [0 1] ...
);

% Load participant info
[basefold, datatype, all_con, ~, ~, participants, ~, ~, activity_tag, ~, ~, ~, ~, ~, ~, ~, ~] = Max_get_param(USING_HPC, get_elec);

% Choose a condition to plot (edit as needed)
condition = char(all_con(1)); % or loop over all_con if you want all conditions

% Call the group plotting function
plot_all_flies_average_2(participants, basefold, datatype, activity_tag, condition, cutoff, timing, plot_cfg);

%Compare CoI metrics (redundancy + synergy) for wake and sleep
% compare_total_coi_across_conditions(participants, basefold, datatype, cutoff, timing)

%Plot similarity between synergy and redundancy patterns
% out_dir = fullfile(basefold, 'Plots', 'Region_Patterns_test')
% compare_region_patterns(participants, basefold, datatype, activity_tag, condition, cutoff, timing, out_dir)