function stats = MI_four_activity_tag_analysis(participants, basefold, ~, activity_tags, cutoff, use_pretrim_data, summary_time_window)
%MI_FOUR_ACTIVITY_TAG_ANALYSIS Compare MI across four activity tags.
%
% Loads single-activity MI files such as:
%   R060721wakeBSLEEP_MI_data.mat
%   R060721sleepBSLEEP_MI_data.mat
%   R060721beginning_sleepBSLEEP_MI_data.mat
%   R060721mid_sleepBSLEEP_MI_data.mat
%
% Then produces:
%   1. Per-electrode MI timecourses with all activity tags.
%   2. Central/peripheral average MI timecourses.
%   3. Average MI boxplots per electrode across activity tags.
%   4. Average MI boxplots for peripheral and central electrodes.

this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir));
addpath(fullfile(project_root, 'Code'));

if nargin < 4 || isempty(activity_tags)
    activity_tags = {'wake', 'sleep', 'active_sleep', 'mid_sleep'};
end
if nargin < 5 || isempty(cutoff)
    cutoff = 10;
end
if nargin < 6 || isempty(use_pretrim_data)
    use_pretrim_data = 0;
end
if nargin < 7
    summary_time_window = [0 74];
end

[default_basefold, ~, all_con, ~, ~, default_participants] = Max_get_comparison_param(0, 0);
if nargin < 1 || isempty(participants)
    participants = default_participants;
end
if nargin < 2 || isempty(basefold)
    basefold = default_basefold;
end
% `datatype` is kept in the signature to mirror related plotting helpers.

times = -25:74;
condition = char(all_con(1));
n_tags = numel(activity_tags);
data_tag = ternary(use_pretrim_data, 'pretrim', 'standard');

comparison_tag = matlab.lang.makeValidName(strjoin(activity_tags, '_'));
save_dir = fullfile(basefold, 'MI_figures', comparison_tag, data_tag, 'four_activity_tags');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

[MI_by_tag, elec_names] = load_four_tag_mi(participants, basefold, condition, activity_tags, use_pretrim_data, times);
n_elec = numel(elec_names);

central_idx = 1:min(cutoff - 1, n_elec);
peripheral_idx = cutoff:n_elec;
region_defs = struct( ...
    'name', {'Peripheral', 'Central'}, ...
    'idx', {peripheral_idx, central_idx});

colors = tag_colors(n_tags);
plot_electrode_timecourses(MI_by_tag, elec_names, times, activity_tags, colors, condition, data_tag, save_dir);
plot_region_timecourses(MI_by_tag, region_defs, elec_names, times, activity_tags, colors, condition, data_tag, save_dir);

summary_mask = make_summary_time_mask(times, summary_time_window);
summary_label = make_summary_window_label(summary_time_window);
electrode_summary = compute_electrode_summary(MI_by_tag, elec_names, summary_mask);
stats.electrodes = plot_electrode_average_boxplot(electrode_summary, elec_names, activity_tags, colors, condition, data_tag, summary_label, save_dir);
stats.regions = plot_region_average_boxplots(electrode_summary, region_defs, elec_names, activity_tags, colors, condition, data_tag, summary_label, save_dir);
end

function [MI_by_tag, elec_names] = load_four_tag_mi(participants, basefold, condition, activity_tags, use_pretrim_data, times)
n_tags = numel(activity_tags);
n_flies = numel(participants);
n_time = numel(times);
MI_by_tag = cell(n_tags, 1);
elec_names = {};

for c = 1:n_tags
    for f = 1:n_flies
        participant = participants{f};
        mi_struct = load_single_activity_mi(basefold, participant, condition, activity_tags{c}, use_pretrim_data);
        if isempty(mi_struct)
            fprintf('Missing MI for %s / %s / %s\n', participant, activity_tags{c}, condition);
            continue
        end

        if isempty(elec_names)
            elec_names = fieldnames(mi_struct);
            for cc = 1:n_tags
                MI_by_tag{cc} = nan(numel(elec_names), n_time, n_flies);
            end
        end

        for e = 1:numel(elec_names)
            elec = elec_names{e};
            if isfield(mi_struct, elec)
                this_mi = mi_struct.(elec);
                MI_by_tag{c}(e, :, f) = this_mi(:)';
            end
        end
    end
end

if isempty(elec_names)
    error('No MI files could be loaded for the requested activity tags.');
end
end

function mi_struct = load_single_activity_mi(basefold, participant, condition, activity_tag, use_pretrim_data)
mi_struct = [];
candidate_names = {sprintf('%s%s%s_MI_data.mat', participant, activity_tag, condition)};
if use_pretrim_data
    candidate_names = [{sprintf('%s_pretrim_%s_%s_MI_data.mat', participant, activity_tag, condition), ...
        sprintf('%s%s%s_MI_data_pretrim.mat', participant, activity_tag, condition)}, candidate_names];
end

for i = 1:numel(candidate_names)
    candidate_path = fullfile(basefold, 'MI_Data', candidate_names{i});
    if ~exist(candidate_path, 'file')
        continue
    end

    S = load(candidate_path, 'MI_stat');
    if isfield(S.MI_stat, participant) && isfield(S.MI_stat.(participant), condition) && isfield(S.MI_stat.(participant).(condition), 'MI')
        mi_struct = S.MI_stat.(participant).(condition).MI;
        return
    end
end
end

function plot_electrode_timecourses(MI_by_tag, elec_names, times, activity_tags, colors, condition, data_tag, save_dir)
n_elec = numel(elec_names);
n_tags = numel(activity_tags);

y_limits = get_global_timecourse_limits(MI_by_tag);
fig = figure('Position', [100 100 1500 650], 'Color', 'w', 'Renderer', 'painters');
tiledlayout(3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

for e = 1:n_elec
    nexttile;
    hold on
    for c = 1:n_tags
        vals = squeeze(MI_by_tag{c}(e, :, :))';
        plot_mean_sem(times, vals, colors(c, :), activity_tags{c});
    end
    xline(0, '--k', 'LineWidth', 0.45);
    title(elec_names{e}, 'Interpreter', 'none');
    xlabel('Time (ms)');
    ylabel('MI (bits)');
    ylim(y_limits);
    grid on
    box off
    if e == 1
        legend(activity_tags, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'none');
    end
    hold off
end

sgtitle(sprintf('MI timecourses across activity tags (%s, %s)', condition, data_tag), 'Interpreter', 'none');
save_figure(fig, save_dir, ['MI_four_activity_tags_timecourses_', condition]);
end

function plot_region_timecourses(MI_by_tag, region_defs, elec_names, times, activity_tags, colors, condition, data_tag, save_dir)
n_tags = numel(activity_tags);
fig = figure('Position', [100 100 1000 430], 'Color', 'w', 'Renderer', 'painters');
tiledlayout(1, numel(region_defs), 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:numel(region_defs)
    nexttile;
    hold on
    region_idx = region_defs(r).idx;
    region_idx = region_idx(region_idx >= 1 & region_idx <= numel(elec_names));
    for c = 1:n_tags
        vals = MI_by_tag{c}(region_idx, :, :);
        vals = reshape(permute(vals, [3 1 2]), [], numel(times));
        plot_mean_sem(times, vals, colors(c, :), activity_tags{c});
    end
    xline(0, '--k', 'LineWidth', 0.45);
    title(region_defs(r).name);
    xlabel('Time (ms)');
    ylabel('MI (bits)');
    grid on
    box off
    if r == 1
        legend(activity_tags, 'Location', 'northeast', 'Box', 'off', 'Interpreter', 'none');
    end
    hold off
end

sgtitle(sprintf('Regional MI timecourses (%s, %s)', condition, data_tag), 'Interpreter', 'none');
save_figure(fig, save_dir, ['MI_four_activity_tags_region_timecourses_', condition]);
end

function summary = compute_electrode_summary(MI_by_tag, elec_names, summary_mask)
n_tags = numel(MI_by_tag);
n_elec = numel(elec_names);
summary = nan(n_elec, n_tags);

for c = 1:n_tags
    for e = 1:n_elec
        vals = MI_by_tag{c}(e, summary_mask, :);
        summary(e, c) = mean(vals(:), 'omitnan');
    end
end
end

function stats_table = plot_electrode_average_boxplot(summary, elec_names, activity_tags, colors, condition, data_tag, summary_label, save_dir)
stats_table = pairwise_summary_stats(summary, activity_tags);
fig = plot_summary_boxplot(summary, activity_tags, colors, ...
    sprintf('Average MI per electrode (%s, %s, %s)', condition, data_tag, summary_label), ...
    'Mean MI per electrode (bits)', stats_table);
save_figure(fig, save_dir, ['MI_four_activity_tags_average_per_electrode_', condition, '_', summary_label]);

summary_table = array2table(summary, 'VariableNames', matlab.lang.makeValidName(activity_tags));
summary_table.Electrode = elec_names(:);
summary_table = movevars(summary_table, 'Electrode', 'Before', 1);
writetable(summary_table, fullfile(save_dir, ['MI_four_activity_tags_average_per_electrode_', condition, '_', summary_label, '.csv']));
writetable(stats_table, fullfile(save_dir, ['MI_four_activity_tags_average_per_electrode_pairwise_stats_', condition, '_', summary_label, '.csv']));
end

function region_stats = plot_region_average_boxplots(summary, region_defs, elec_names, activity_tags, colors, condition, data_tag, summary_label, save_dir)
region_stats = struct();
for r = 1:numel(region_defs)
    region_idx = region_defs(r).idx;
    region_idx = region_idx(region_idx >= 1 & region_idx <= numel(elec_names));
    region_summary = summary(region_idx, :);
    stats_table = pairwise_summary_stats(region_summary, activity_tags);
    region_stats.(matlab.lang.makeValidName(region_defs(r).name)) = stats_table;

    fig = plot_summary_boxplot(region_summary, activity_tags, colors, ...
        sprintf('%s electrode average MI (%s, %s, %s)', region_defs(r).name, condition, data_tag, summary_label), ...
        sprintf('Mean MI per %s electrode (bits)', lower(region_defs(r).name)), stats_table);
    save_figure(fig, save_dir, ['MI_four_activity_tags_average_', matlab.lang.makeValidName(region_defs(r).name), '_', condition, '_', summary_label]);
    writetable(stats_table, fullfile(save_dir, ['MI_four_activity_tags_average_', matlab.lang.makeValidName(region_defs(r).name), '_pairwise_stats_', condition, '_', summary_label, '.csv']));
end
end

function fig = plot_summary_boxplot(summary, activity_tags, colors, plot_title, y_label, stats_table)
n_tags = numel(activity_tags);
fig = figure('Color', 'w', 'Position', [100 100 680 540], 'Renderer', 'painters');
hold on
all_vals = summary(isfinite(summary));

for c = 1:n_tags
    vals = summary(:, c);
    vals = vals(isfinite(vals));
    boxchart(c * ones(size(vals)), vals, ...
        'BoxFaceColor', colors(c, :), ...
        'BoxFaceAlpha', 0.58, ...
        'MarkerStyle', 'none', ...
        'BoxWidth', 0.40, ...
        'LineWidth', 0.40);
    scatter(c + 0.045*randn(size(vals)), vals, 35, colors(c, :), ...
        'filled', 'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'none');
end

xlim([0.45, n_tags + 0.55]);
xticks(1:n_tags);
xticklabels(activity_tags);
xtickangle(25);
ylabel(y_label);
title(plot_title, 'Interpreter', 'none');
grid on
box off
set(gca, 'FontSize', 12, 'LineWidth', 0.40);

if ~isempty(all_vals)
    add_pairwise_brackets(stats_table, activity_tags, all_vals);
end
hold off
end

function stats_table = pairwise_summary_stats(summary, activity_tags)
pairs = nchoosek(1:numel(activity_tags), 2);
stats_table = table('Size', [size(pairs, 1) 9], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double','logical'}, ...
    'VariableNames', {'ActivityA','ActivityB','N','MeanDiff_BMinusA','PairedT','PairedDf','PairedP','PairedPFDR','FDR_sig'});

for p = 1:size(pairs, 1)
    c1 = pairs(p, 1);
    c2 = pairs(p, 2);
    vals_a = summary(:, c1);
    vals_b = summary(:, c2);
    valid = isfinite(vals_a) & isfinite(vals_b);
    diff_vals = vals_b(valid) - vals_a(valid);

    stats_table.ActivityA(p) = string(activity_tags{c1});
    stats_table.ActivityB(p) = string(activity_tags{c2});
    stats_table.N(p) = sum(valid);
    stats_table.MeanDiff_BMinusA(p) = mean(diff_vals, 'omitnan');

    if sum(valid) >= 2
        try
            [~, p_val, ~, t_stats] = ttest(vals_a(valid), vals_b(valid));
            stats_table.PairedT(p) = t_stats.tstat;
            stats_table.PairedDf(p) = t_stats.df;
            stats_table.PairedP(p) = p_val;
        catch
            stats_table.PairedT(p) = NaN;
            stats_table.PairedDf(p) = NaN;
            stats_table.PairedP(p) = NaN;
        end
    end
end

[stats_table.PairedPFDR, stats_table.FDR_sig] = benjamini_hochberg(stats_table.PairedP, 0.05);
end

function add_pairwise_brackets(stats_table, activity_tags, all_vals)
y_min = min(all_vals);
y_max = max(all_vals);
y_range = max(y_max - y_min, eps);
ylim([y_min - 0.06*y_range, y_max + 0.45*y_range]);

sig_rows = find(stats_table.FDR_sig & isfinite(stats_table.PairedPFDR));
if isempty(sig_rows)
    text(mean(xlim), y_max + 0.10*y_range, 'no FDR-significant pairwise differences', ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.35 0.35 0.35]);
    return
end

for i = 1:min(numel(sig_rows), 6)
    row = stats_table(sig_rows(i), :);
    x1 = find(strcmp(activity_tags, char(row.ActivityA)), 1);
    x2 = find(strcmp(activity_tags, char(row.ActivityB)), 1);
    y = y_max + (0.07 + 0.055*i) * y_range;
    plot([x1 x1 x2 x2], [y y+0.014*y_range y+0.014*y_range y], 'k-', 'LineWidth', 0.40);
    text(mean([x1 x2]), y + 0.018*y_range, sprintf('q=%.2g', row.PairedPFDR), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end
end

function plot_mean_sem(times, vals, color, display_name)
times = times(:)';
if isvector(vals)
    vals = vals(:)';
end
mean_vals = mean(vals, 1, 'omitnan');
n_vals = sum(isfinite(vals), 1);
sem_vals = std(vals, 0, 1, 'omitnan') ./ sqrt(n_vals);
sem_vals(n_vals == 0) = NaN;
mean_vals = mean_vals(:)';
sem_vals = sem_vals(:)';
fill([times, fliplr(times)], [mean_vals + sem_vals, fliplr(mean_vals - sem_vals)], ...
    color, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(times, mean_vals, 'Color', color, 'LineWidth', 0.75, 'DisplayName', display_name);
end

function y_limits = get_global_timecourse_limits(MI_by_tag)
all_vals = [];
for c = 1:numel(MI_by_tag)
    all_vals = [all_vals; MI_by_tag{c}(:)]; %#ok<AGROW>
end
all_vals = all_vals(isfinite(all_vals));
if isempty(all_vals)
    y_limits = [0 1];
else
    y_range = max(range(all_vals), eps);
    y_limits = [min(all_vals) - 0.05*y_range, max(all_vals) + 0.05*y_range];
end
end

function colors = tag_colors(n_tags)
base_colors = [0.20 0.45 0.80; 0.85 0.33 0.10; 0.35 0.65 0.35; 0.55 0.35 0.75; 0.65 0.55 0.20];
colors = base_colors(1:n_tags, :);
end

function summary_mask = make_summary_time_mask(times, summary_time_window)
if isempty(summary_time_window)
    summary_mask = true(size(times));
else
    summary_mask = times >= summary_time_window(1) & times <= summary_time_window(2);
end
end

function label = make_summary_window_label(summary_time_window)
if isempty(summary_time_window)
    label = 'full_epoch';
else
    label = sprintf('%dto%dms', summary_time_window(1), summary_time_window(2));
end
end

function save_figure(fig, save_dir, filename_stem)
saveas(fig, fullfile(save_dir, [filename_stem, '.fig']));
saveas(fig, fullfile(save_dir, [filename_stem, '.png']));
print_editable_svg(fig, fullfile(save_dir, [filename_stem, '.svg']));
close(fig);
end

function print_editable_svg(fig_handle, svg_file)
thin_figure_for_illustrator(fig_handle);
set(fig_handle, 'Renderer', 'painters');
print(fig_handle, svg_file, '-dsvg', '-vector');
end

function thin_figure_for_illustrator(fig_handle)
thinLineWidth = 0.40;
objects = findall(fig_handle);
for objInd = 1:numel(objects)
    obj = objects(objInd);
    if isprop(obj, 'LineWidth')
        try
            currentWidth = get(obj, 'LineWidth');
            if isnumeric(currentWidth) && isscalar(currentWidth) && currentWidth > thinLineWidth
                set(obj, 'LineWidth', thinLineWidth);
            end
        catch
        end
    end
end

axesObjects = findall(fig_handle, 'Type', 'axes');
for axInd = 1:numel(axesObjects)
    try
        set(axesObjects(axInd), 'LineWidth', thinLineWidth);
    catch
    end
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function [adj_p, sig] = benjamini_hochberg(pvals, alpha)
adj_p = nan(size(pvals));
sig = false(size(pvals));
valid = ~isnan(pvals);
if ~any(valid)
    return
end

pv = pvals(valid);
[sorted_p, sort_idx] = sort(pv);
rank = (1:numel(sorted_p))';
adj_sorted = sorted_p .* numel(sorted_p) ./ rank;
for i = numel(adj_sorted)-1:-1:1
    adj_sorted(i) = min(adj_sorted(i), adj_sorted(i+1));
end
adj_sorted(adj_sorted > 1) = 1;
tmp = nan(size(pv));
tmp(sort_idx) = adj_sorted;
adj_p(valid) = tmp;
sig(valid) = adj_p(valid) <= alpha;
end
