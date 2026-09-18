function stats = compare_total_coi_four_activity_tags(participants, basefold, datatype, cutoff, times, activity_tags)
%COMPARE_TOTAL_COI_FOUR_ACTIVITY_TAGS Pairwise total CoI comparisons across activity tags.
%
% This is a four-condition companion to compare_total_coi_across_conditions.
% It uses the same total-CoI extraction logic: for each fly, region, and
% activity tag, it loads CoI data, keeps significant redundancy or synergy
% values using mask_r/mask_s, then pools those values across electrode pairs.
%
% Defaults compare:
%   wake, sleep, beginning_sleep, mid_sleep
%
% Outputs are saved under:
%   <basefold>/Results/Group_Analysis/<activity_tag_1>_<...>_<activity_tag_n>/

if nargin < 6 || isempty(activity_tags)
    activity_tags = {'wake', 'sleep', 'active_sleep', 'mid_sleep'};
end

if numel(activity_tags) < 2
    error('Provide at least two activity tags.');
end

region_fields = {'peripheral', 'central', 'central_peripheral', 'peripheral_peripheral', 'central_central'};
region_labels = {'Peripheral', 'Central', 'Central-Peripheral', 'Peripheral-Peripheral', 'Central-Central'};
activity_tag = make_activity_comparison_tag(activity_tags);
save_dir = fullfile(basefold, 'Results', 'Group_Analysis', activity_tag);
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

feature_defs = struct( ...
    'key', {'redundancy', 'synergy'}, ...
    'mask_field', {'mask_r', 'mask_s'}, ...
    'name', {'Redundancy', 'Synergy'} ...
);

stats = struct();
for k = 1:numel(feature_defs)
    stats.(feature_defs(k).key) = run_feature_analysis( ...
        feature_defs(k).name, feature_defs(k).mask_field, participants, basefold, datatype, ...
        cutoff, times, activity_tags, region_fields, region_labels, save_dir);
end
end

function stats_table = run_feature_analysis(feature_name, mask_field, participants, basefold, datatype, ...
    cutoff, times, activity_tags, region_fields, region_labels, save_dir)

n_regions = numel(region_fields);
n_tags = numel(activity_tags);
n_flies = numel(participants);

[~, ~, all_con, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = Max_get_param(0, 0);
condition = char(all_con(1));

feature_all = cell(n_regions, n_tags);
fly_means = cell(n_regions, n_tags);

for c = 1:n_tags
    this_activity_tag = activity_tags{c};
    for r = 1:n_regions
        region_name = region_fields{r};
        all_vals = [];
        means_per_fly = nan(n_flies, 1);

        for f = 1:n_flies
            participant_name = char(participants(f));
            CoIData = get_cached_plotting_CoI_raw_rf(basefold, datatype, participant_name, this_activity_tag, condition, cutoff, times);

            if ~isfield(CoIData, region_name) || ~isfield(CoIData.(region_name), 'FFi_all') || ~isfield(CoIData.(region_name), mask_field)
                continue
            end

            mask_vals = CoIData.(region_name).(mask_field);
            FFi = CoIData.(region_name).FFi_all;
            vals = FFi(mask_vals > 0);
            vals = vals(:);
            vals = vals(isfinite(vals));

            all_vals = [all_vals; vals]; %#ok<AGROW>
            if ~isempty(vals)
                means_per_fly(f) = mean(vals, 'omitnan');
            end
        end

        feature_all{r, c} = all_vals;
        fly_means{r, c} = means_per_fly;
    end
end

pairs = nchoosek(1:n_tags, 2);
n_rows = n_regions * size(pairs, 1);
stats_table = table('Size', [n_rows 19], ...
    'VariableTypes', {'string','string','string','string','double','double','double','double','double','double','double','double','double','double','double','double','double','double','logical'}, ...
    'VariableNames', {'Region','Label','ActivityA','ActivityB','N_A_AllPairs','N_B_AllPairs','MeanA_AllPairs','MeanB_AllPairs','Diff_BMinusA_AllPairs','WelchT','WelchDf','WelchP','RanksumP','CohensD_AllPairs','N_Flies_Paired','Diff_BMinusA_PerFly','PairedT','PairedP','FDR_sig'});

row_ind = 0;
for r = 1:n_regions
    for p = 1:size(pairs, 1)
        c1 = pairs(p, 1);
        c2 = pairs(p, 2);
        row_ind = row_ind + 1;

        vals_a = clean_values(feature_all{r, c1});
        vals_b = clean_values(feature_all{r, c2});
        means_a = fly_means{r, c1};
        means_b = fly_means{r, c2};
        valid_fly = isfinite(means_a) & isfinite(means_b);
        diff_fly = means_b(valid_fly) - means_a(valid_fly);

        [welch_t, welch_df, welch_p, ranksum_p, cohens_d] = compare_unpaired_values(vals_a, vals_b);
        [paired_t, paired_p] = compare_paired_fly_means(means_a(valid_fly), means_b(valid_fly));

        stats_table.Region(row_ind) = string(region_fields{r});
        stats_table.Label(row_ind) = string(region_labels{r});
        stats_table.ActivityA(row_ind) = string(activity_tags{c1});
        stats_table.ActivityB(row_ind) = string(activity_tags{c2});
        stats_table.N_A_AllPairs(row_ind) = numel(vals_a);
        stats_table.N_B_AllPairs(row_ind) = numel(vals_b);
        stats_table.MeanA_AllPairs(row_ind) = mean(vals_a, 'omitnan');
        stats_table.MeanB_AllPairs(row_ind) = mean(vals_b, 'omitnan');
        stats_table.Diff_BMinusA_AllPairs(row_ind) = mean(vals_b, 'omitnan') - mean(vals_a, 'omitnan');
        stats_table.WelchT(row_ind) = welch_t;
        stats_table.WelchDf(row_ind) = welch_df;
        stats_table.WelchP(row_ind) = welch_p;
        stats_table.RanksumP(row_ind) = ranksum_p;
        stats_table.CohensD_AllPairs(row_ind) = cohens_d;
        stats_table.N_Flies_Paired(row_ind) = sum(valid_fly);
        stats_table.Diff_BMinusA_PerFly(row_ind) = mean(diff_fly, 'omitnan');
        stats_table.PairedT(row_ind) = paired_t;
        stats_table.PairedP(row_ind) = paired_p;

        fprintf('%s %s: %s vs %s, Welch t(%g)=%.3f, p=%.4g, ranksum p=%.4g, d=%.3f\n', ...
            region_labels{r}, upper(feature_name), activity_tags{c1}, activity_tags{c2}, ...
            welch_df, welch_t, welch_p, ranksum_p, cohens_d);
    end
end

[stats_table.WelchPFDR, stats_table.FDR_sig] = benjamini_hochberg(stats_table.WelchP, 0.05);
writetable(stats_table, fullfile(save_dir, [feature_name '_FourActivityTag_PairwiseStats.csv']));

for r = 1:n_regions
    region_rows = strcmp(stats_table.Region, region_fields{r});
    plot_four_tag_region_boxplot(feature_all(r, :), feature_name, activity_tags, region_labels{r}, ...
        stats_table(region_rows, :), save_dir);
end

plot_compact_allpairs_summary(feature_all, feature_name, activity_tags, region_fields, region_labels, stats_table, save_dir);
end

function vals = clean_values(vals)
vals = vals(:);
vals = vals(isfinite(vals));
end

function [t_stat, t_df, p_val, p_ranksum, cohens_d] = compare_unpaired_values(vals_a, vals_b)
if isempty(vals_a) || isempty(vals_b)
    t_stat = NaN;
    t_df = NaN;
    p_val = NaN;
    p_ranksum = NaN;
    cohens_d = NaN;
    return
end

[~, p_val, ~, t_stats] = ttest2(vals_a, vals_b, 'Vartype', 'unequal');
t_stat = t_stats.tstat;
t_df = t_stats.df;

try
    p_ranksum = ranksum(vals_a, vals_b);
catch
    p_ranksum = NaN;
end

pooled_sd = sqrt(((numel(vals_a)-1)*var(vals_a,'omitnan') + (numel(vals_b)-1)*var(vals_b,'omitnan')) / ...
    max(numel(vals_a)+numel(vals_b)-2, 1));
if pooled_sd == 0 || isnan(pooled_sd)
    cohens_d = NaN;
else
    cohens_d = (mean(vals_a,'omitnan') - mean(vals_b,'omitnan')) / pooled_sd;
end
end

function [t_stat, p_val] = compare_paired_fly_means(vals_a, vals_b)
if numel(vals_a) < 2 || numel(vals_b) < 2
    t_stat = NaN;
    p_val = NaN;
    return
end

try
    [~, p_val, ~, t_stats] = ttest(vals_a, vals_b);
    t_stat = t_stats.tstat;
catch
    t_stat = NaN;
    p_val = NaN;
end
end

function plot_four_tag_region_boxplot(region_feature_all, feature_name, activity_tags, region_label, region_stats, save_dir)
colors = [0.20 0.45 0.80; 0.85 0.33 0.10; 0.35 0.65 0.35; 0.55 0.35 0.75; 0.65 0.55 0.20];
max_scatter = 900;
n_tags = numel(activity_tags);
activity_labels = activity_display_labels(activity_tags);

fig = figure('Color', 'w', 'Position', [100 100 980 560]);
hold on

all_plot_vals = [];
for c = 1:n_tags
    vals = clean_values(region_feature_all{c});
    if isempty(vals)
        continue
    end

    color_ind = 1 + mod(c - 1, size(colors, 1));
    boxchart(c * ones(size(vals)), vals, ...
        'BoxFaceColor', colors(color_ind, :), ...
        'BoxFaceAlpha', 0.58, ...
        'MarkerStyle', 'none', ...
        'BoxWidth', 0.42, ...
        'LineWidth', 0.40);

    if numel(vals) > max_scatter
        vals_for_scatter = vals(randperm(numel(vals), max_scatter));
    else
        vals_for_scatter = vals;
    end
    scatter(c + 0.055*randn(size(vals_for_scatter)), vals_for_scatter, ...
        8, colors(color_ind, :), 'filled', 'MarkerFaceAlpha', 0.14, ...
        'MarkerEdgeColor', 'none');
    all_plot_vals = [all_plot_vals; vals]; %#ok<AGROW>
end

if isempty(all_plot_vals)
    close(fig);
    warning('No %s values available for %s.', feature_name, region_label);
    return
end

xlim([0.45, n_tags + 0.55]);
xticks(1:n_tags);
xticklabels(activity_labels);
xtickangle(25);
ylabel([feature_name, ' (bits)']);
title(sprintf('%s: %s across conditions', region_label, feature_name), 'Interpreter', 'none');
grid on
box off
set(gca, 'FontSize', 11, 'LineWidth', 0.40);

add_significance_brackets(region_stats, activity_tags, all_plot_vals, feature_name);

comparison_tag = make_activity_comparison_tag(activity_tags);
region_tag = matlab.lang.makeValidName(region_label);
saveas(fig, fullfile(save_dir, [feature_name '_FourActivityTags_' comparison_tag '_' region_tag '.png']));
saveas(fig, fullfile(save_dir, [feature_name '_FourActivityTags_' comparison_tag '_' region_tag '.fig']));
print_editable_svg(fig, fullfile(save_dir, [feature_name '_FourActivityTags_' comparison_tag '_' region_tag '.svg']));
close(fig);
end

function add_significance_brackets(region_stats, activity_tags, all_plot_vals, feature_name)
[~, y_max, y_range] = set_zoomed_ylim(all_plot_vals, 1.18, percentile_limits_for_feature(feature_name));

sig_rows = find(region_stats.FDR_sig & isfinite(region_stats.WelchPFDR));
if isempty(sig_rows)
    return
end

max_brackets = min(numel(sig_rows), 6);
for i = 1:max_brackets
    row = region_stats(sig_rows(i), :);
    x1 = find(strcmp(activity_tags, char(row.ActivityA)), 1);
    x2 = find(strcmp(activity_tags, char(row.ActivityB)), 1);
    if isempty(x1) || isempty(x2)
        continue
    end

    y = y_max + (0.08 + 0.095*i) * y_range;
    plot([x1 x1 x2 x2], [y y+0.015*y_range y+0.015*y_range y], ...
        'k-', 'LineWidth', 0.55, 'Clipping', 'off');
    text(mean([x1 x2]), y + 0.006*y_range, p_to_star(row.WelchPFDR), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'Clipping', 'off');
end
end

function plot_compact_allpairs_summary(feature_all, feature_name, activity_tags, region_fields, region_labels, stats_table, save_dir)
% Compact all-pairs summary matching compare_total_coi_across_conditions.
%
% Only Peripheral, Central, and Central-Peripheral are shown. Within each
% region, the boxes contain all significant redundancy/synergy values pooled
% across electrode pairs and flies.

region_inds_to_plot = 1:min(3, numel(region_labels));
n_regions_to_plot = numel(region_inds_to_plot);
n_tags = numel(activity_tags);
activity_labels = activity_display_labels(activity_tags);
colors = [0.20 0.45 0.80; 0.85 0.33 0.10; 0.35 0.65 0.35; 0.55 0.35 0.75; 0.65 0.55 0.20];
offsets = linspace(-0.30, 0.30, n_tags);
max_scatter = 500;

fig = figure('Color', 'w', 'Position', [100 100 1300 560]);
hold on

all_plot_vals = [];
for plot_r = 1:n_regions_to_plot
    r = region_inds_to_plot(plot_r);
    for c = 1:n_tags
        vals = clean_values(feature_all{r, c});
        if isempty(vals)
            continue
        end

        color_ind = 1 + mod(c - 1, size(colors, 1));
        x = plot_r + offsets(c);
        boxchart(x * ones(size(vals)), vals, ...
            'BoxFaceColor', colors(color_ind, :), ...
            'BoxFaceAlpha', 0.58, ...
            'MarkerStyle', 'none', ...
            'BoxWidth', 0.12, ...
            'LineWidth', 0.40);

        if numel(vals) > max_scatter
            vals_for_scatter = vals(randperm(numel(vals), max_scatter));
        else
            vals_for_scatter = vals;
        end
        scatter(x + 0.020*randn(size(vals_for_scatter)), vals_for_scatter, ...
            6, colors(color_ind, :), 'filled', 'MarkerFaceAlpha', 0.12, ...
            'MarkerEdgeColor', 'none');

        all_plot_vals = [all_plot_vals; vals]; %#ok<AGROW>
    end
end

if isempty(all_plot_vals)
    close(fig);
    warning('No all-pairs %s values available for compact summary.', feature_name);
    return
end

[~, y_max, y_range] = set_zoomed_ylim(all_plot_vals, 1.42, percentile_limits_for_feature(feature_name));

for plot_r = 1:n_regions_to_plot
    r = region_inds_to_plot(plot_r);
    region_rows = strcmp(stats_table.Region, region_fields{r});
    region_stats = stats_table(region_rows, :);
    add_compact_region_brackets(region_stats, activity_tags, offsets, plot_r, y_max, y_range);
end

xlim([0.45, n_regions_to_plot + 0.55]);
xticks(1:n_regions_to_plot);
xticklabels(region_labels(region_inds_to_plot));
xtickangle(20);
ylabel([feature_name, ' (bits)']);
title([feature_name, ' across conditions'], 'Interpreter', 'none', 'FontWeight', 'bold');
grid on
box off
set(gca, 'FontSize', 11, 'LineWidth', 0.40);

legend_handles = gobjects(n_tags, 1);
for c = 1:n_tags
    color_ind = 1 + mod(c - 1, size(colors, 1));
    legend_handles(c) = scatter(nan, nan, 30, colors(color_ind, :), 'filled', 'MarkerEdgeColor', 'none');
end
legend(legend_handles, activity_labels, 'Location', 'northeastoutside', 'Box', 'off', 'Interpreter', 'none');

comparison_tag = make_activity_comparison_tag(activity_tags);
saveas(fig, fullfile(save_dir, [feature_name '_FourActivityTags_AllPairs_CompactSummary_' comparison_tag '.png']));
saveas(fig, fullfile(save_dir, [feature_name '_FourActivityTags_AllPairs_CompactSummary_' comparison_tag '.fig']));
print_editable_svg(fig, fullfile(save_dir, [feature_name '_FourActivityTags_AllPairs_CompactSummary_' comparison_tag '.svg']));
close(fig);
end

function add_compact_region_brackets(region_stats, activity_tags, offsets, region_x, y_max, y_range)
sig_rows = find(region_stats.FDR_sig & isfinite(region_stats.WelchPFDR));
if isempty(sig_rows)
    return
end

max_brackets = min(numel(sig_rows), 6);
for i = 1:max_brackets
    row = region_stats(sig_rows(i), :);
    c1 = find(strcmp(activity_tags, char(row.ActivityA)), 1);
    c2 = find(strcmp(activity_tags, char(row.ActivityB)), 1);
    if isempty(c1) || isempty(c2)
        continue
    end

    x1 = region_x + offsets(c1);
    x2 = region_x + offsets(c2);
    y = y_max + (0.07 + 0.090*i) * y_range;
    plot([x1 x1 x2 x2], [y y+0.014*y_range y+0.014*y_range y], ...
        'k-', 'LineWidth', 0.40, 'Clipping', 'off');
    text(mean([x1 x2]), y + 0.005*y_range, p_to_star(row.WelchPFDR), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Clipping', 'off');
end
end

function [y_min, y_max, y_range] = set_zoomed_ylim(all_plot_vals, top_padding_fraction, percentile_limits)
% Use feature-specific robust limits so rare CoI values do not flatten the
% paper figure. Redundancy mainly needs upper-tail clipping; synergy is
% negative, so it mainly needs lower-tail clipping.
if nargin < 3 || isempty(percentile_limits)
    percentile_limits = [0 99];
end
finite_vals = all_plot_vals(isfinite(all_plot_vals));
if isempty(finite_vals)
    y_min = 0;
    y_max = 1;
    y_range = 1;
    ylim([0 1]);
    return
end

if numel(finite_vals) >= 20
    if percentile_limits(1) <= 0
        y_min = min(finite_vals);
    else
        y_min = prctile(finite_vals, percentile_limits(1));
    end
    if percentile_limits(2) >= 100
        y_max = max(finite_vals);
    else
        y_max = prctile(finite_vals, percentile_limits(2));
    end
else
    y_min = min(finite_vals);
    y_max = max(finite_vals);
end

if y_min == y_max
    y_range = max(abs(y_max), 1);
    y_min = y_min - 0.5*y_range;
    y_max = y_max + 0.5*y_range;
else
    y_range = y_max - y_min;
end

ylim([y_min - 0.06*y_range, y_max + top_padding_fraction*y_range]);
end

function percentile_limits = percentile_limits_for_feature(feature_name)
if strcmpi(feature_name, 'Synergy')
    percentile_limits = [0.5 99];
else
    percentile_limits = [0 99];
end
end

function star_txt = p_to_star(p_val)
if ~isfinite(p_val)
    star_txt = '';
elseif p_val < 0.001
    star_txt = '***';
elseif p_val < 0.01
    star_txt = '**';
elseif p_val < 0.05
    star_txt = '*';
else
    star_txt = '';
end
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

function labels = activity_display_labels(activity_tags)
labels = activity_tags;
for label_ind = 1:numel(labels)
    switch labels{label_ind}
        case 'mid_sleep'
            labels{label_ind} = 'deep sleep';
        case 'active_sleep'
            labels{label_ind} = 'active sleep';
        otherwise
            labels{label_ind} = strrep(labels{label_ind}, '_', ' ');
    end
end
end
function comparison_tag = make_activity_comparison_tag(activity_tags)
comparison_tag = matlab.lang.makeValidName(strjoin(activity_tags, '_'));
end

function CoIData = get_cached_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, times)
persistent CoIDataCache
if isempty(CoIDataCache)
    CoIDataCache = containers.Map('KeyType','char','ValueType','any');
end

cutoff_tag = strrep(mat2str(cutoff), ' ', '_');
cache_key = strjoin({ ...
    basefold, datatype, participant_name, activity_tag, condition, cutoff_tag, ...
    ['nTimes' num2str(numel(times))]}, '|');

if isKey(CoIDataCache, cache_key)
    CoIData = CoIDataCache(cache_key);
    return
end

CoIData = max_get_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, times);
CoIDataCache(cache_key) = CoIData;
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
