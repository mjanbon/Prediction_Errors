%% MI ANALYSIS: Activity State Comparison with Permutation-Based Significance
% This script loads mutual information (MI) comparison data computed by
% main_MI_ERP_comparison.m and visualizes the differences between two
% activity states (e.g., wake vs sleep, beginning_sleep vs mid_sleep, etc.)
% across all electrodes.
%
% The significance testing uses cluster-based permutation tests (implemented
% in cnm_MI_stimtime_perm_diff.m) rather than traditional parametric tests.
%
% Author: Maxime Janbon
% Date: February 2026

%% SETUP AND PARAMETERS
% Load MI comparison data with significance masks from permutation testing
USING_HPC = 0;  % 0 = local machine, 1 = Cambridge HPC, 2 = QMUL HPC
get_elec = 0;   % 0 = don't get electrodes of interest (not needed for plotting)
use_pretrim_data = 1; % 0 = load standard files, 1 = load *_pretrim files
times = -25:74; % Time vector in ms (baseline: -25 to 0, post-stim: 0 to 74)
stim_idx = 26;  % Index of stimulus onset (time = 0 ms)
n_flies = 8;    % Total number of flies in the dataset

data_tags = {'standard', 'pretrim'};
data_tag = data_tags{1 + logical(use_pretrim_data)};
summary_time_window = [0 74]; % ms; use [] to average across the full plotted epoch

%% LOAD PARAMETERS
% Get experimental parameters using the comparison parameter function
% This loads parameters for BOTH activity states (activity_tag_1 and activity_tag_2)
% Examples: wake vs sleep, beginning_sleep vs mid_sleep, etc.
[basefold, datatype, all_con, condition,subject, participants, EoI,...
srate, activity_tag_1, deviant_group_number_1, standard_group_number_1,...
activity_tag_2, deviant_group_number_2, standard_group_number_2, corrected,...
stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_comparison_param(USING_HPC, get_elec);

% Create comparison string for filenames and titles
comparison_str = strcat(activity_tag_1, '_', activity_tag_2);
save_dir = fullfile(basefold, 'MI_figures', comparison_str, data_tag);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% INITIALIZE DATA STRUCTURES
% Preallocate structures to store MI data for all participants
% Structure format: struct.(participant).(condition).(electrode) = time_vector
MI_diff_data = struct();  % MI difference (activity_tag_1 - activity_tag_2)
MI_1_data = struct();     % MI for activity_tag_1 (first activity state)
MI_2_data = struct();     % MI for activity_tag_2 (second activity state)
sigMask_data = struct();  % Significance mask from permutation testing (binary: 0 or ±1)

%% LOAD MI COMPARISON DATA
% Load comparison data for all participants and conditions
% Data files are created by main_MI_ERP_comparison.m which computes:
%   - MI_1: Mutual information during activity_tag_1 state
%   - MI_2: Mutual information during activity_tag_2 state  
%   - MI_diff: Difference between activity_tag_1 and activity_tag_2 MI
%   - sigMask: Significance mask from cluster-based permutation test
for i = 1:length(participants)
    for con = 1:length(all_con)
        participant_id = participants{i};
        con_name = all_con{con};

        MI_name_candidates = {sprintf('%s%s%s_MI_data.mat', participant_id, comparison_str, con_name)};
        if use_pretrim_data
            % Preferred pretrim pattern:
            % 'R290721_pretrim_wake_sleep_BSLEEP_MI_data.mat'
            MI_name_candidates = [{
                sprintf('%s_pretrim_%s_%s_MI_data.mat', participant_id, comparison_str, con_name), ...
                sprintf('%s%s%s_MI_data_pretrim.mat', participant_id, comparison_str, con_name)}, ...
                MI_name_candidates];
        end

        MI_name = '';
        MI_path = '';
        for name_idx = 1:length(MI_name_candidates)
            candidate_path = fullfile(basefold, 'MI_Data', MI_name_candidates{name_idx});
            if exist(candidate_path, 'file')
                MI_name = MI_name_candidates{name_idx};
                MI_path = candidate_path;
                break;
            end
        end

        if isempty(MI_name)
            fprintf('Warning: File not found. Tried: %s\n', strjoin(MI_name_candidates, ' | '));
            continue;
        end
        
        % Load MI statistics from saved file
        load(MI_path, 'MI_stat');
        
        % Extract channel names (e.g., E1, E2, ..., E15)
        chans = fieldnames(MI_stat.MI_diff);
        
        % Store data for each electrode
        for ch = 1:length(chans)
            chan_name = chans{ch};
            MI_diff_data.(participants{i}).(all_con{con}).(chan_name) = MI_stat.MI_diff.(chan_name);
            MI_1_data.(participants{i}).(all_con{con}).(chan_name) = MI_stat.MI_1.(chan_name);
            MI_2_data.(participants{i}).(all_con{con}).(chan_name) = MI_stat.MI_2.(chan_name);
            sigMask_data.(participants{i}).(all_con{con}).(chan_name) = MI_stat.sigMask.(chan_name);
        end
    end
end

%% GET ELECTRODE INFORMATION
% Extract electrode names from the first participant's data after loading
% Assumes all participants have the same electrode configuration
sub_condition = all_con{1};  % e.g., 'BSLEEP'
elec_names = fieldnames(MI_1_data.(participants{1}).(sub_condition));
n_elec = length(elec_names);  % Typically 15 electrodes

%% PREPARE DATA MATRICES FOR PLOTTING
% Organize MI data into matrices for statistical analysis and visualization
% Matrix dimensions: [n_flies × n_time_points]

% Get electrode names
sub_condition = all_con{1};
elec_names = fieldnames(MI_1_data.(participants{1}).(sub_condition));
n_elec = length(elec_names);

% --- Prepare data for plotting: MI_1 (wake) and MI_2 (sleep) ---
n_time = length(times);
MI_1_mat_all = cell(n_elec, 1);
MI_2_mat_all = cell(n_elec, 1);
sigMask_all = cell(n_elec, 1);

% Calculate global y-limits
ymin = inf; ymax = -inf;

for ch = 1:n_elec
    chan_name = elec_names{ch};
    MI_1_mat = nan(n_flies, n_time);
    MI_2_mat = nan(n_flies, n_time);
    
    for i = 1:n_flies
        if isfield(MI_1_data.(participants{i}).(sub_condition), chan_name)
            mi_1 = MI_1_data.(participants{i}).(sub_condition).(chan_name);
            mi_2 = MI_2_data.(participants{i}).(sub_condition).(chan_name);
            MI_1_mat(i,:) = mi_1(:)';
            MI_2_mat(i,:) = mi_2(:)';
        end
    end
    
    % Calculate means and SEMs
    mean_MI_1 = mean(MI_1_mat, 1, 'omitnan');
    mean_MI_2 = mean(MI_2_mat, 1, 'omitnan');
    n_nonan_1 = sum(~isnan(MI_1_mat), 1);
    n_nonan_2 = sum(~isnan(MI_2_mat), 1);
    sem_MI_1 = std(MI_1_mat, 0, 1, 'omitnan') ./ sqrt(n_nonan_1);
    sem_MI_2 = std(MI_2_mat, 0, 1, 'omitnan') ./ sqrt(n_nonan_2);
    
    % Update global min/max
    ymin = min([ymin, min(mean_MI_1 - sem_MI_1), min(mean_MI_2 - sem_MI_2)]);
    ymax = max([ymax, max(mean_MI_1 + sem_MI_1), max(mean_MI_2 + sem_MI_2)]);
    
    % Store matrices for later plotting
    MI_1_mat_all{ch} = MI_1_mat;
    MI_2_mat_all{ch} = MI_2_mat;
    
    % Load significance masks for this electrode across all flies
    sigMask_mat = nan(n_flies, n_time);
    for i = 1:n_flies
        if isfield(sigMask_data.(participants{i}).(sub_condition), chan_name)
            sig = sigMask_data.(participants{i}).(sub_condition).(chan_name);
            sigMask_mat(i,:) = sig(:)';
        end
    end
    
    % Create consensus significance mask:
    % A time point is significant if >50% of flies show significance
    sigMask_all{ch} = mean(sigMask_mat, 1, 'omitnan') >= 0.5;
end

%% PLOT: PAPER-STYLE SUMMARY ACROSS ELECTRODES
% Each point is one electrode. For each electrode, MI is averaged across
% flies and across summary_time_window. The two boxplots compare the two
% activity tags directly, with paired electrode lines.

if isempty(summary_time_window)
    summary_time_mask = true(size(times));
    summary_window_label = 'full_epoch';
    summary_window_title = 'full epoch';
else
    summary_time_mask = times >= summary_time_window(1) & times <= summary_time_window(2);
    summary_window_label = sprintf('%dto%dms', summary_time_window(1), summary_time_window(2));
    summary_window_title = sprintf('%d to %d ms', summary_time_window(1), summary_time_window(2));
end

summary_MI_1 = nan(n_elec,1);
summary_MI_2 = nan(n_elec,1);
for ch = 1:n_elec
    summary_MI_1(ch) = mean(MI_1_mat_all{ch}(:,summary_time_mask), 'all', 'omitnan');
    summary_MI_2(ch) = mean(MI_2_mat_all{ch}(:,summary_time_mask), 'all', 'omitnan');
end

valid_summary = isfinite(summary_MI_1) & isfinite(summary_MI_2);
summary_MI_1_plot = summary_MI_1(valid_summary);
summary_MI_2_plot = summary_MI_2(valid_summary);
summary_elec_names = elec_names(valid_summary);

[~, summary_p, ~, summary_stats] = ttest(summary_MI_1_plot, summary_MI_2_plot);
summary_diff = summary_MI_2_plot - summary_MI_1_plot;

summary_fig = figure('Color','w','Position',[100 100 520 520]);
hold on
boxchart(ones(size(summary_MI_1_plot)), summary_MI_1_plot, ...
    'BoxFaceColor', [0 0.4470 0.7410], 'MarkerStyle', 'none');
boxchart(2*ones(size(summary_MI_2_plot)), summary_MI_2_plot, ...
    'BoxFaceColor', [0.8500 0.3250 0.0980], 'MarkerStyle', 'none');
jitter_width = 0.04;
summary_x_1 = 1 + jitter_width*randn(size(summary_MI_1_plot));
summary_x_2 = 2 + jitter_width*randn(size(summary_MI_2_plot));
for elec_idx = 1:numel(summary_MI_1_plot)
    plot([summary_x_1(elec_idx), summary_x_2(elec_idx)], ...
        [summary_MI_1_plot(elec_idx), summary_MI_2_plot(elec_idx)], ...
        '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8);
end
scatter(summary_x_1, summary_MI_1_plot, 36, ...
    [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.65);
scatter(summary_x_2, summary_MI_2_plot, 36, ...
    [0.8500 0.3250 0.0980], 'filled', 'MarkerFaceAlpha', 0.65);
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({activity_tag_1, activity_tag_2});
ylabel('Mean MI across electrodes (bits)');
title({['Average MI per electrode: ', sub_condition, ', ', data_tag], ...
    ['Window: ', summary_window_title, ...
    ', paired t(', num2str(summary_stats.df), ') = ', num2str(summary_stats.tstat,'%.2f'), ...
    ', p = ', num2str(summary_p,'%.3g')]}, ...
    'Interpreter','none');
grid on
box off
set(gca,'FontSize',12,'LineWidth',1);
hold off

summary_table_var_names = matlab.lang.makeValidName({'Electrode', activity_tag_1, activity_tag_2, [activity_tag_2 '_minus_' activity_tag_1]});
summary_table = table(summary_elec_names(:), summary_MI_1_plot(:), summary_MI_2_plot(:), summary_diff(:), ...
    'VariableNames', summary_table_var_names);
writetable(summary_table, fullfile(save_dir, ['MI_summary_per_electrode_', sub_condition, '_', summary_window_label, '.csv']));
saveas(summary_fig, fullfile(save_dir, ['MI_summary_per_electrode_', sub_condition, '_', summary_window_label, '.fig']));
saveas(summary_fig, fullfile(save_dir, ['MI_summary_per_electrode_', sub_condition, '_', summary_window_label, '.png']));
print(summary_fig, fullfile(save_dir, ['MI_summary_per_electrode_', sub_condition, '_', summary_window_label, '.svg']), '-dsvg');

%% DISPLAY SIGNIFICANCE RESULTS
% Print summary of significant time points for each electrode
% based on cluster-based permutation testing
%
% Significance determined by:
%   1. Cluster-based permutation test in cnm_MI_stimtime_perm_diff
%   2. Consensus threshold: >=50% of flies show significance at that time point

fprintf('\n=== Significance Masks from Permutation Testing ===\n');
fprintf('Data source: %s\n', data_tag);
fprintf('Threshold: p < 0.05 (cluster-corrected)\n');
fprintf('Consensus: Time points where >50%% of flies show significance\n\n');

for ch = 1:n_elec
    chan_name = elec_names{ch};
    sig_consensus = sigMask_all{ch};
    sig_times_idx = find(sig_consensus);
    
    fprintf('Electrode: %s\n', chan_name);
    if ~isempty(sig_times_idx)
        fprintf('  Significant time points: ');
        fprintf('%d ', times(sig_times_idx));
        fprintf('ms\n');
        fprintf('  Duration: %d ms\n', length(sig_times_idx));
    else
        fprintf('  No significant differences found.\n');
    end
end

%% PLOT: MI COMPARISON WITH SIGNIFICANCE OVERLAY
% Tiled layout showing MI time courses for both activity states with
% black dots indicating significant time points where activity_tag_1 and
% activity_tag_2 MI significantly differ (cluster-corrected p < 0.05)
%
% Black dots are positioned near the top of each subplot to maximize visibility
% Dots indicate consensus significance (>50% of flies)

figure('Position',[100 100 1500 600]);
tiledlayout(3,5,'TileSpacing','compact','Padding','compact');
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980]; % activity 1 (blue), activity 2 (orange)

for ch = 1:n_elec
    chan_name = elec_names{ch};
    
    % Retrieve data and significance mask
    MI_1_mat = MI_1_mat_all{ch};
    MI_2_mat = MI_2_mat_all{ch};
    sig_consensus = sigMask_all{ch};
    
    % Calculate mean and SEM
    mean_MI_1 = mean(MI_1_mat, 1, 'omitnan');
    mean_MI_2 = mean(MI_2_mat, 1, 'omitnan');
    n_nonan_1 = sum(~isnan(MI_1_mat), 1);
    n_nonan_2 = sum(~isnan(MI_2_mat), 1);
    sem_MI_1 = std(MI_1_mat, 0, 1, 'omitnan') ./ sqrt(n_nonan_1);
    sem_MI_2 = std(MI_2_mat, 0, 1, 'omitnan') ./ sqrt(n_nonan_2);
    
    nexttile;
    
    % Plot activity_tag_1 MI with SEM shading
    fill([times, fliplr(times)], ...
         [mean_MI_1 + sem_MI_1, fliplr(mean_MI_1 - sem_MI_1)], ...
         colors(1,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
        h1 = plot(times, mean_MI_1, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', activity_tag_1);
    
    % Plot activity_tag_2 MI with SEM shading
    fill([times, fliplr(times)], ...
         [mean_MI_2 + sem_MI_2, fliplr(mean_MI_2 - sem_MI_2)], ...
         colors(2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        h2 = plot(times, mean_MI_2, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', activity_tag_2);
    
    % Overlay significance markers as black dots near the top
    sig_times = times(sig_consensus);
    if ~isempty(sig_times)
        y_top = ymax - 0.05 * (ymax - ymin);  % Position dots 5% below top
        plot(sig_times, y_top * ones(size(sig_times)), 'k.', 'MarkerSize', 12);
    end
    
    % Add stimulus onset marker
    xline(0, '--k', 'LineWidth', 1);
    
    % Labels and formatting
    title(chan_name);
    xlabel('Time (ms)');
    ylabel('MI (bits)');
    ylim([ymin ymax]);
    
    % Show legend only on first subplot
    if ch == 1
        legend([h1 h2], {activity_tag_1, activity_tag_2}, 'Location', 'northeast', 'Box', 'off');
    end
end

sgtitle(['MI comparison with significance (', sub_condition, ', ', data_tag, ') - black dots: permutation test p<0.05']);

% Save the figure in multiple formats
saveas(gcf, fullfile(save_dir, ['MI_timecourse_comparison_with_significance_', sub_condition, '.fig']));
saveas(gcf, fullfile(save_dir, ['MI_timecourse_comparison_with_significance_', sub_condition, '.png']));
print(gcf, fullfile(save_dir, ['MI_timecourse_comparison_with_significance_', sub_condition, '.svg']), '-dsvg');

%% PLOT: MI DIFFERENCE WITH SIGNIFICANCE OVERLAY
% Show the difference between activity states (activity_tag_2 - activity_tag_1) for all electrodes
% with black dots indicating significant time points from cluster-based permutation test
%
% Positive values: higher MI in activity_tag_2
% Negative values: higher MI in activity_tag_1

figure('Position',[100 100 1500 600]);
tiledlayout(3,5,'TileSpacing','compact','Padding','compact');

% Calculate global y-limits for MI difference
ymin_diff = inf; ymax_diff = -inf;
for ch = 1:n_elec
    chan_name = elec_names{ch};
    MI_diff_mat = nan(n_flies, n_time);
    
    for i = 1:n_flies
        if isfield(MI_diff_data.(participants{i}).(sub_condition), chan_name)
            mi_diff = MI_diff_data.(participants{i}).(sub_condition).(chan_name);
            MI_diff_mat(i,:) = mi_diff(:)';
        end
    end
    
    mean_MI_diff = mean(MI_diff_mat, 1, 'omitnan');
    n_nonan = sum(~isnan(MI_diff_mat), 1);
    sem_MI_diff = std(MI_diff_mat, 0, 1, 'omitnan') ./ sqrt(n_nonan);
    ymin_diff = min([ymin_diff, min(mean_MI_diff - sem_MI_diff)]);
    ymax_diff = max([ymax_diff, max(mean_MI_diff + sem_MI_diff)]);
end

for ch = 1:n_elec
    chan_name = elec_names{ch};
    
    % Load MI difference from saved data
    MI_diff_mat = nan(n_flies, n_time);
    for i = 1:n_flies
        if isfield(MI_diff_data.(participants{i}).(sub_condition), chan_name)
            mi_diff = MI_diff_data.(participants{i}).(sub_condition).(chan_name);
            MI_diff_mat(i,:) = mi_diff(:)';
        end
    end
    
    sig_consensus = sigMask_all{ch};
    
    % Calculate mean and SEM of difference
    mean_MI_diff = mean(MI_diff_mat, 1, 'omitnan');
    n_nonan = sum(~isnan(MI_diff_mat), 1);
    sem_MI_diff = std(MI_diff_mat, 0, 1, 'omitnan') ./ sqrt(n_nonan);
    
    nexttile;
    
    % Plot MI difference with SEM shading
    fill([times, fliplr(times)], ...
         [mean_MI_diff + sem_MI_diff, fliplr(mean_MI_diff - sem_MI_diff)], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    plot(times, mean_MI_diff, 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
    
    % Add zero line
    yline(0, '--k', 'LineWidth', 0.5);
    
    % Overlay significance markers as black dots near the top
    sig_times = times(sig_consensus);
    if ~isempty(sig_times)
        y_top = ymax_diff - 0.05 * (ymax_diff - ymin_diff);
        plot(sig_times, y_top * ones(size(sig_times)), 'k.', 'MarkerSize', 12);
    end
    
    % Add stimulus onset marker
    xline(0, '--r', 'LineWidth', 1);
    
    % Labels and formatting
    title(chan_name);
    xlabel('Time (ms)');
    ylabel([activity_tag_2, ' - ', activity_tag_1, ' (bits)']);
    ylim([ymin_diff ymax_diff]);
end

sgtitle(['MI difference: ', activity_tag_2, ' - ', activity_tag_1, ' (', sub_condition, ', ', data_tag, ') - black dots: p<0.05']);

% Save the figure
saveas(gcf, fullfile(save_dir, ['MI_difference_with_significance_', sub_condition, '.fig']));
saveas(gcf, fullfile(save_dir, ['MI_difference_with_significance_', sub_condition, '.png']));
print(gcf, fullfile(save_dir, ['MI_difference_with_significance_', sub_condition, '.svg']), '-dsvg');
