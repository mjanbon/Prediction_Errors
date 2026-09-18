function stats = compare_total_coi_across_conditions(participants, basefold, datatype, cutoff, times, activity_tags)
%COMPARE_TOTAL_COI_ACROSS_CONDITIONS
% Run the SAME plotting/stats pipeline for both redundancy and synergy.
%
% Keeping one implementation ensures figure/style edits affect both analyses.

if nargin < 6 || isempty(activity_tags)
    activity_tags = {'wake', 'sleep'};
end

region_fields = {'peripheral', 'central', 'central_peripheral', 'peripheral_peripheral', 'central_central'};
region_labels = {'Peripheral', 'Central', 'Central-Peripheral', 'Peripheral-Peripheral', 'Central-Central'};
activity_comparison_tag = make_activity_comparison_tag(activity_tags);
save_dir = fullfile(basefold, 'Results', 'Group_Analysis', activity_comparison_tag);
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

feature_defs = struct( ...
    'key', {'redundancy', 'synergy'}, ...
    'mask_field', {'mask_r', 'mask_s'}, ...
    'name', {'Redundancy', 'Synergy'} ...
);

stats = struct();
compact_summary_data = struct();
for k = 1:numel(feature_defs)
    [stats_per_fly, stats_all_pairs, feature_all] = run_feature_analysis(feature_defs(k).key, feature_defs(k).name, feature_defs(k).mask_field, ...
        participants, basefold, datatype, cutoff, times, activity_tags, region_fields, region_labels, save_dir);
    stats.(feature_defs(k).key).per_fly = stats_per_fly;
    stats.(feature_defs(k).key).all_pairs = stats_all_pairs;
    compact_summary_data.(feature_defs(k).key).name = feature_defs(k).name;
    compact_summary_data.(feature_defs(k).key).feature_all = feature_all;
    compact_summary_data.(feature_defs(k).key).stats_allpairs = stats_all_pairs;
end

if isfield(compact_summary_data, 'redundancy') && isfield(compact_summary_data, 'synergy')
    plot_combined_compact_distribution_summary(compact_summary_data, activity_tags, region_labels, save_dir);
end
end

function [stats, stats_allpairs, feature_all] = run_feature_analysis(feature_key, feature_name, mask_field, ...
    participants, basefold, datatype, cutoff, times, activity_tags, region_fields, region_labels, save_dir)

n_regions = numel(region_fields);
n_conditions = numel(activity_tags);
n_flies = numel(participants);
activity_labels = activity_display_labels(activity_tags);

% Store all values for each region/condition
feature_all = cell(n_regions, n_conditions);
fly_means = cell(n_regions, n_conditions);

[~, ~, all_con, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = Max_get_param(0, 0);
condition = char(all_con(1));

for c = 1:n_conditions
    activity_tag = activity_tags{c};
    for r = 1:n_regions
        region_name = region_fields{r};
        all_vals = [];
        means_per_fly = nan(n_flies, 1);
        for f = 1:n_flies
            participant_name = char(participants(f));
            CoIData = get_cached_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, times);
            if isfield(CoIData, region_name) && isfield(CoIData.(region_name), 'FFi_all') && isfield(CoIData.(region_name), mask_field)
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
        end
        feature_all{r, c} = all_vals;
        fly_means{r, c} = means_per_fly;
    end
end

stats = table('Size', [n_regions 10], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Region','Label','N_Flies','WakeMean','SleepMean','DiffMean','tStat','tDf','pValue','CohensD'});

stats_allpairs = table('Size', [n_regions 18], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Region','Label','N_WakePairs','N_SleepPairs','WakeMean','SleepMean','DiffMean','WakeMedian','SleepMedian','DiffMedian','WakeAbsMedian','SleepAbsMedian','DiffAbsMedian','tStat','tDf','pValue','RanksumP','CohensD'});

for r = 1:n_regions
    region_label = region_labels{r};
    all_wake = feature_all{r,1};
    all_sleep = feature_all{r,2};
    all_wake = all_wake(isfinite(all_wake));
    all_sleep = all_sleep(isfinite(all_sleep));

    % ---------- All-pairs analysis ----------
    if ~isempty(all_wake) && ~isempty(all_sleep)
        [~, p_allpairs, ~, stats_all] = ttest2(all_wake, all_sleep, 'Vartype', 'unequal');
        t_stat_all = stats_all.tstat;
        t_df_all = stats_all.df;
        try
            p_ranksum = ranksum(all_wake, all_sleep);
        catch
            p_ranksum = NaN;
        end

        pooled_sd = sqrt(((numel(all_wake)-1)*var(all_wake,'omitnan') + (numel(all_sleep)-1)*var(all_sleep,'omitnan')) / ...
            max(numel(all_wake)+numel(all_sleep)-2, 1));
        if pooled_sd == 0 || isnan(pooled_sd)
            cohens_d_all = NaN;
        else
            cohens_d_all = (mean(all_wake,'omitnan') - mean(all_sleep,'omitnan')) / pooled_sd;
        end

        stats_allpairs.Region(r) = string(region_fields{r});
        stats_allpairs.Label(r) = string(region_label);
        stats_allpairs.N_WakePairs(r) = numel(all_wake);
        stats_allpairs.N_SleepPairs(r) = numel(all_sleep);
        stats_allpairs.WakeMean(r) = mean(all_wake, 'omitnan');
        stats_allpairs.SleepMean(r) = mean(all_sleep, 'omitnan');
        stats_allpairs.DiffMean(r) = mean(all_sleep, 'omitnan') - mean(all_wake, 'omitnan');
        stats_allpairs.WakeMedian(r) = median(all_wake, 'omitnan');
        stats_allpairs.SleepMedian(r) = median(all_sleep, 'omitnan');
        stats_allpairs.DiffMedian(r) = median(all_sleep, 'omitnan') - median(all_wake, 'omitnan');
        stats_allpairs.WakeAbsMedian(r) = median(abs(all_wake), 'omitnan');
        stats_allpairs.SleepAbsMedian(r) = median(abs(all_sleep), 'omitnan');
        stats_allpairs.DiffAbsMedian(r) = median(abs(all_sleep), 'omitnan') - median(abs(all_wake), 'omitnan');
        stats_allpairs.tStat(r) = t_stat_all;
        stats_allpairs.tDf(r) = t_df_all;
        stats_allpairs.pValue(r) = p_allpairs;
        stats_allpairs.RanksumP(r) = p_ranksum;
        stats_allpairs.CohensD(r) = cohens_d_all;

        fprintf('%s %s (all pairs): t(%d)=%.3f, p=%.4g, ranksum p=%.4g, Cohen''s d=%.3f, Nwake=%d, Nsleep=%d\n', ...
            region_label, upper(feature_key), t_df_all, t_stat_all, p_allpairs, p_ranksum, cohens_d_all, numel(all_wake), numel(all_sleep));

        n_boot = 1000;
        diff_boot = nan(n_boot, 1);
        n_w = numel(all_wake);
        n_s = numel(all_sleep);
        for b = 1:n_boot
            w_b = all_wake(randi(n_w, n_w, 1));
            s_b = all_sleep(randi(n_s, n_s, 1));
            diff_boot(b) = mean(s_b, 'omitnan') - mean(w_b, 'omitnan');
        end

        fig_all = figure('Color', 'w', 'Position', [100 100 1000 450], 'Renderer', 'painters');
        tl_all = tiledlayout(fig_all, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl_all, sprintf('%s: %s vs %s %s', region_label, activity_labels{1}, activity_labels{2}, feature_name), 'FontWeight', 'bold', 'Interpreter', 'none');

        nexttile;
        hold on;
        boxchart(ones(size(all_wake)), all_wake, 'BoxFaceColor', [0.20 0.45 0.80], 'MarkerStyle', 'none', 'LineWidth', 0.45);
        boxchart(2*ones(size(all_sleep)), all_sleep, 'BoxFaceColor', [0.85 0.33 0.10], 'MarkerStyle', 'none', 'LineWidth', 0.45);
        max_scatter = 1500;
        if numel(all_wake) > max_scatter
            wake_plot = all_wake(randperm(numel(all_wake), max_scatter));
        else
            wake_plot = all_wake;
        end
        if numel(all_sleep) > max_scatter
            sleep_plot = all_sleep(randperm(numel(all_sleep), max_scatter));
        else
            sleep_plot = all_sleep;
        end
        scatter(1 + 0.05*randn(size(wake_plot)), wake_plot, 10, [0.20 0.45 0.80], 'filled', 'MarkerFaceAlpha', 0.20, 'MarkerEdgeColor', 'none');
        scatter(2 + 0.05*randn(size(sleep_plot)), sleep_plot, 10, [0.85 0.33 0.10], 'filled', 'MarkerFaceAlpha', 0.20, 'MarkerEdgeColor', 'none');
        xlim([0.5 2.5]);
        xticks([1 2]);
        xticklabels(activity_labels);
        ylabel(sprintf('%s (bits)', feature_name));
        grid on; box off;
        set(gca, 'FontSize', 11, 'LineWidth', 0.45);
        y_pair = [all_wake; all_sleep];
        y_top = max(y_pair) + 0.03 * max(range(y_pair), eps);
        text(1.5, y_top, sprintf('Welch t: t(%d)=%.2f, p=%.3g', t_df_all, t_stat_all, p_allpairs), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
        hold off;

        nexttile;
        hold on;
        boxchart(ones(size(diff_boot)), diff_boot, 'BoxFaceColor', [0.35 0.35 0.35], 'MarkerStyle', 'none', 'LineWidth', 0.45);
        scatter(1 + 0.05*randn(size(diff_boot)), diff_boot, 8, [0.1 0.1 0.1], 'filled', 'MarkerFaceAlpha', 0.10, 'MarkerEdgeColor', 'none');
        yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
        xlim([0.5 1.5]);
        xticks(1);
        xticklabels({[activity_labels{2}, '-', activity_labels{1}, ' (boot mean)']});
        ylabel('Mean difference distribution');
        grid on; box off;
        set(gca, 'FontSize', 11, 'LineWidth', 0.45);
        ci_boot = prctile(diff_boot, [2.5 97.5]);
        text(1, max(diff_boot) + 0.05*max(range(diff_boot), eps), ...
            sprintf('mean=%.3f, 95%% CI=[%.3f, %.3f]\nd=%.3f', ...
            mean(diff_boot,'omitnan'), ci_boot(1), ci_boot(2), cohens_d_all), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
        hold off;

        saveas(fig_all, fullfile(save_dir, [feature_name '_AllPairs_' make_activity_comparison_tag(activity_tags) '_' region_label '.png']));
        saveas(fig_all, fullfile(save_dir, [feature_name '_AllPairs_' make_activity_comparison_tag(activity_tags) '_' region_label '.fig']));
        print_editable_svg(fig_all, fullfile(save_dir, [feature_name '_AllPairs_' make_activity_comparison_tag(activity_tags) '_' region_label '.svg']));
        close(fig_all);
    else
        stats_allpairs.Region(r) = string(region_fields{r});
        stats_allpairs.Label(r) = string(region_label);
        stats_allpairs.N_WakePairs(r) = numel(all_wake);
        stats_allpairs.N_SleepPairs(r) = numel(all_sleep);
        stats_allpairs.WakeMean(r) = NaN;
        stats_allpairs.SleepMean(r) = NaN;
        stats_allpairs.DiffMean(r) = NaN;
        stats_allpairs.WakeMedian(r) = NaN;
        stats_allpairs.SleepMedian(r) = NaN;
        stats_allpairs.DiffMedian(r) = NaN;
        stats_allpairs.WakeAbsMedian(r) = NaN;
        stats_allpairs.SleepAbsMedian(r) = NaN;
        stats_allpairs.DiffAbsMedian(r) = NaN;
        stats_allpairs.tStat(r) = NaN;
        stats_allpairs.tDf(r) = NaN;
        stats_allpairs.pValue(r) = NaN;
        stats_allpairs.RanksumP(r) = NaN;
        stats_allpairs.CohensD(r) = NaN;
        fprintf('%s %s (all pairs): insufficient data for statistical comparison.\n', region_label, upper(feature_key));
    end

    % ---------- Per-fly paired analysis ----------
    fly_means_wake = fly_means{r,1};
    fly_means_sleep = fly_means{r,2};
    diff_per_fly = fly_means_sleep - fly_means_wake;
    valid_diff = isfinite(diff_per_fly);
    diff_per_fly = diff_per_fly(valid_diff);
    fly_means_wake = fly_means_wake(valid_diff);
    fly_means_sleep = fly_means_sleep(valid_diff);

    [~, p_paired, ~, paired_stats] = ttest(fly_means_wake, fly_means_sleep);
    t_stat = paired_stats.tstat;
    t_df = paired_stats.df;

    if isempty(diff_per_fly) || std(diff_per_fly) == 0
        cohens_d = NaN;
    else
        cohens_d = mean(diff_per_fly, 'omitnan') / std(diff_per_fly, 'omitnan');
    end

    try
        p_signedrank = signrank(fly_means_wake, fly_means_sleep);
    catch
        p_signedrank = NaN;
    end

    stats.Region(r) = string(region_fields{r});
    stats.Label(r) = string(region_label);
    stats.N_Flies(r) = numel(diff_per_fly);
    stats.WakeMean(r) = mean(fly_means_wake, 'omitnan');
    stats.SleepMean(r) = mean(fly_means_sleep, 'omitnan');
    stats.DiffMean(r) = mean(diff_per_fly, 'omitnan');
    stats.tStat(r) = t_stat;
    stats.tDf(r) = t_df;
    stats.pValue(r) = p_paired;
    stats.CohensD(r) = cohens_d;

    fprintf('%s %s (per-fly paired): t(%d)=%.3f, p=%.4g, signrank p=%.4g, Cohen''s d=%.3f\n', ...
        region_label, upper(feature_key), t_df, t_stat, p_paired, p_signedrank, cohens_d);

    fig = figure('Color', 'w', 'Position', [100 100 900 450], 'Renderer', 'painters');
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: %s - %s %s', region_label, activity_labels{1}, activity_labels{2}, feature_name), 'FontWeight', 'bold', 'Interpreter', 'none');

    nexttile;
    hold on;
    boxchart(ones(size(fly_means_wake)), fly_means_wake, 'BoxFaceColor', [0.20 0.45 0.80], 'MarkerStyle', 'none', 'LineWidth', 0.45);
    boxchart(2*ones(size(fly_means_sleep)), fly_means_sleep, 'BoxFaceColor', [0.85 0.33 0.10], 'MarkerStyle', 'none', 'LineWidth', 0.45);
    scatter(1 + 0.06*randn(size(fly_means_wake)), fly_means_wake, 36, [0.20 0.45 0.80], 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
    scatter(2 + 0.06*randn(size(fly_means_sleep)), fly_means_sleep, 36, [0.85 0.33 0.10], 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels(activity_labels);
    ylabel(sprintf('Per-fly mean %s', lower(feature_name)));
    grid on; box off;
    set(gca, 'FontSize', 11, 'LineWidth', 0.45);
    text(1.5, max([fly_means_wake; fly_means_sleep]) + 0.02*range([fly_means_wake; fly_means_sleep]), ...
        sprintf('paired t: t(%d)=%.2f, p=%.3g', t_df, t_stat, p_paired), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
    hold off;

    nexttile;
    hold on;
    boxchart(ones(size(diff_per_fly)), diff_per_fly, 'BoxFaceColor', [0.40 0.40 0.40], 'MarkerStyle', 'none', 'LineWidth', 0.45);
    scatter(1 + 0.06*randn(size(diff_per_fly)), diff_per_fly, 42, [0.15 0.15 0.15], 'filled', 'MarkerFaceAlpha', 0.70, 'MarkerEdgeColor', 'none');
    yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    xlim([0.5 1.5]);
    xticks(1);
    xticklabels({[activity_labels{2}, '-', activity_labels{1}]});
    ylabel(sprintf('Difference in per-fly mean %s', lower(feature_name)));
    grid on; box off;
    set(gca, 'FontSize', 11, 'LineWidth', 0.45);
    text(1, max(diff_per_fly) + 0.05*range(diff_per_fly + eps), ...
        sprintf('mean diff = %.3f\nCohen''s d = %.3f', mean(diff_per_fly, 'omitnan'), cohens_d), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
    hold off;

    saveas(fig, fullfile(save_dir, [feature_name '_' make_activity_difference_tag(activity_tags) '_' region_label '.png']));
    saveas(fig, fullfile(save_dir, [feature_name '_' make_activity_difference_tag(activity_tags) '_' region_label '.fig']));
    print_editable_svg(fig, fullfile(save_dir, [feature_name '_' make_activity_difference_tag(activity_tags) '_' region_label '.svg']));
    close(fig);
end
% ---------- Multiple comparisons correction (FDR) ----------
alpha_fdr = 0.05;
% prepare columns
stats_allpairs.pFDR = nan(n_regions,1);
stats_allpairs.FDR_sig = false(n_regions,1);
stats.pFDR = nan(n_regions,1);
stats.FDR_sig = false(n_regions,1);

[adj_all, sig_all] = benjamini_hochberg(stats_allpairs.pValue, alpha_fdr);
stats_allpairs.pFDR = adj_all;
stats_allpairs.FDR_sig = sig_all;

[adj_pf, sig_pf] = benjamini_hochberg(stats.pValue, alpha_fdr);
stats.pFDR = adj_pf;
stats.FDR_sig = sig_pf;

fprintf('Applied Benjamini-Hochberg FDR (alpha=%.2f): %d/%d all-pairs significant, %d/%d per-fly significant\n', ...
    alpha_fdr, sum(sig_all,'omitnan'), n_regions, sum(sig_pf,'omitnan'), n_regions);

plot_compact_allpairs_summary(feature_all, feature_name, activity_tags, region_labels, stats_allpairs, save_dir);

writetable(stats, fullfile(save_dir, [feature_name '_' make_activity_difference_tag(activity_tags) '_Stats.csv']));
writetable(stats_allpairs, fullfile(save_dir, [feature_name '_AllPairs_' make_activity_comparison_tag(activity_tags) '_Stats.csv']));
fprintf('Saved %s per-fly stats to %s\n', lower(feature_name), fullfile(save_dir, [feature_name '_' make_activity_difference_tag(activity_tags) '_Stats.csv']));
fprintf('Saved %s all-pairs stats to %s\n', lower(feature_name), fullfile(save_dir, [feature_name '_AllPairs_' make_activity_comparison_tag(activity_tags) '_Stats.csv']));
end

function plot_compact_allpairs_summary(feature_all, feature_name, activity_tags, region_labels, stats_allpairs, save_dir)
% Compact paper-style all-pairs summary.
%
% Each box contains all significant redundancy/synergy values from all
% electrode pairs and all flies for one activity state and one region.

region_inds_to_plot = 1:min(3,numel(region_labels)); % Peripheral, Central, Central-Peripheral
region_labels_to_plot = region_labels(region_inds_to_plot);
n_regions_to_plot = numel(region_inds_to_plot);
activity_labels = activity_display_labels(activity_tags);
colors = [0.20 0.45 0.80; 0.85 0.33 0.10];
offsets = [-0.18, 0.18];
max_scatter = 600;
n_boot = 1000;

fig = figure('Color', 'w', 'Position', [100 100 1250 430], 'Renderer', 'painters');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, [feature_name, ': ', activity_labels{1}, ' vs ', activity_labels{2}, ' summary'], 'Interpreter', 'none');

all_plot_vals = [];
diff_mean = nan(n_regions_to_plot,1);
diff_ci = nan(n_regions_to_plot,2);

nexttile(1)
hold on
for plot_r = 1:n_regions_to_plot
    r = region_inds_to_plot(plot_r);
    for c = 1:2
        vals = feature_all{r,c};
        vals = vals(:);
        vals = vals(isfinite(vals));
        x = plot_r + offsets(c);

        if isempty(vals)
            continue
        end

        boxchart(x * ones(size(vals)), vals, ...
            'BoxFaceColor', colors(c,:), ...
            'BoxFaceAlpha', 0.55, ...
            'MarkerStyle', 'none', ...
            'BoxWidth', 0.26, ...
            'LineWidth', 0.45);

        if numel(vals) > max_scatter
            vals_for_scatter = vals(randperm(numel(vals), max_scatter));
        else
            vals_for_scatter = vals;
        end

        scatter(x + 0.035*randn(size(vals_for_scatter)), vals_for_scatter, ...
            8, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.14, ...
            'MarkerEdgeColor', 'none');

        all_plot_vals = [all_plot_vals; vals]; %#ok<AGROW>
    end

    all_wake = feature_all{r,1};
    all_sleep = feature_all{r,2};
    all_wake = all_wake(isfinite(all_wake));
    all_sleep = all_sleep(isfinite(all_sleep));
    if ~isempty(all_wake) && ~isempty(all_sleep)
        diff_mean(plot_r) = mean(all_sleep,'omitnan') - mean(all_wake,'omitnan');
        diff_boot = nan(n_boot,1);
        for boot_ind = 1:n_boot
            wake_boot = all_wake(randi(numel(all_wake),numel(all_wake),1));
            sleep_boot = all_sleep(randi(numel(all_sleep),numel(all_sleep),1));
            diff_boot(boot_ind) = mean(sleep_boot,'omitnan') - mean(wake_boot,'omitnan');
        end
        diff_ci(plot_r,:) = prctile(diff_boot,[2.5 97.5]);
    end
end

if isempty(all_plot_vals)
    close(fig);
    warning('No all-pairs %s values available for compact summary.', feature_name);
    return
end

[y_min, y_max, y_range] = set_zoomed_ylim(all_plot_vals, 2.7);
y_limits = [y_min, y_max];
y_pad = 0.08 * max(y_range, eps);

for plot_r = 1:n_regions_to_plot
    r = region_inds_to_plot(plot_r);
    y_text = y_limits(2) + 0.38*y_pad;
    stats_txt = sprintf('Welch''s t(%d)=%.2f\np=%.2g, d=%.2f', ...
        stats_allpairs.tDf(r), ...
        stats_allpairs.tStat(r), ...
        stats_allpairs.pValue(r), ...
        stats_allpairs.CohensD(r));
    text(plot_r, y_text, stats_txt, 'HorizontalAlignment', 'center', 'FontSize', 8);

    if is_significant_allpairs(stats_allpairs, r)
        x1 = plot_r + offsets(1);
        x2 = plot_r + offsets(2);
        y_bracket = y_limits(2) + 1.35*y_pad;
        y_tick = 0.12*y_pad;
        plot([x1 x1 x2 x2], ...
            [y_bracket-y_tick, y_bracket, y_bracket, y_bracket-y_tick], ...
            'k-', 'LineWidth', 0.5);
        text(plot_r, y_bracket + 0.08*y_pad, p_to_star(stats_allpairs.pValue(r)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 13, ...
            'FontWeight', 'bold');
    end
end

xlim([0.45, n_regions_to_plot + 0.55]);
xticks(1:n_regions_to_plot);
xticklabels(region_labels_to_plot);
xtickangle(25);
ylabel([feature_name, ' (bits)']);
title('Condition distributions');
grid on
box off
set(gca, 'FontSize', 11, 'LineWidth', 0.45);

legend_handles = gobjects(2,1);
for c = 1:2
    legend_handles(c) = scatter(nan, nan, 30, colors(c,:), 'filled', 'MarkerEdgeColor', 'none');
end
legend(legend_handles, activity_labels(1:2), 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'Box', 'off');

nexttile(2)
hold on
yline(0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.5);
for plot_r = 1:n_regions_to_plot
    if isnan(diff_mean(plot_r))
        continue
    end
    err_low = diff_mean(plot_r) - diff_ci(plot_r,1);
    err_high = diff_ci(plot_r,2) - diff_mean(plot_r);
    errorbar(plot_r, diff_mean(plot_r), err_low, err_high, ...
        'o', ...
        'Color', [0.2 0.2 0.2], ...
        'MarkerFaceColor', [0.2 0.2 0.2], ...
        'MarkerEdgeColor', [0.2 0.2 0.2], ...
        'LineWidth', 0.65, ...
        'CapSize', 9);
end
xlim([0.45, n_regions_to_plot + 0.55]);
xticks(1:n_regions_to_plot);
xticklabels(region_labels_to_plot);
xtickangle(25);
ylabel([activity_labels{2}, ' - ', activity_labels{1}, ' (bits)']);
title('Mean Difference');
grid on
box off
set(gca, 'FontSize', 11, 'LineWidth', 0.45);
hold off

comparison_tag = make_activity_comparison_tag(activity_tags);
saveas(fig, fullfile(save_dir, [feature_name '_AllPairs_CompactSummary_' comparison_tag '.png']));
saveas(fig, fullfile(save_dir, [feature_name '_AllPairs_CompactSummary_' comparison_tag '.fig']));
print_editable_svg(fig, fullfile(save_dir, [feature_name '_AllPairs_CompactSummary_' comparison_tag '.svg']));
close(fig);
end

function plot_combined_compact_distribution_summary(compact_summary_data, activity_tags, region_labels, save_dir)
% Paper summary with only the condition-distribution panels: redundancy on
% the left, synergy on the right.
comparison_tag = make_activity_comparison_tag(activity_tags);
fig = figure('Color', 'w', 'Position', [100 100 1250 430], 'Renderer', 'painters');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'CoI across conditions', 'Interpreter', 'none', 'FontWeight', 'bold');

axRed = nexttile(1);
draw_compact_condition_distribution_panel(axRed, ...
    compact_summary_data.redundancy.feature_all, ...
    compact_summary_data.redundancy.name, ...
    activity_tags, region_labels, ...
    compact_summary_data.redundancy.stats_allpairs, true);

axSyn = nexttile(2);
draw_compact_condition_distribution_panel(axSyn, ...
    compact_summary_data.synergy.feature_all, ...
    compact_summary_data.synergy.name, ...
    activity_tags, region_labels, ...
    compact_summary_data.synergy.stats_allpairs, false);

saveas(fig, fullfile(save_dir, ['CoI_AllPairs_CompactDistributionSummary_' comparison_tag '.png']));
saveas(fig, fullfile(save_dir, ['CoI_AllPairs_CompactDistributionSummary_' comparison_tag '.fig']));
print_editable_svg(fig, fullfile(save_dir, ['CoI_AllPairs_CompactDistributionSummary_' comparison_tag '.svg']));
close(fig);
end

function draw_compact_condition_distribution_panel(ax, feature_all, feature_name, activity_tags, region_labels, stats_allpairs, show_legend)
region_inds_to_plot = 1:min(3,numel(region_labels)); % Peripheral, Central, Central-Peripheral
region_labels_to_plot = region_labels(region_inds_to_plot);
n_regions_to_plot = numel(region_inds_to_plot);
activity_labels = activity_display_labels(activity_tags);
colors = [0.20 0.45 0.80; 0.85 0.33 0.10];
offsets = [-0.18, 0.18];
max_scatter = 600;

axes(ax);
hold(ax, 'on');
all_plot_vals = [];
for plot_r = 1:n_regions_to_plot
    r = region_inds_to_plot(plot_r);
    for c = 1:2
        vals = feature_all{r,c};
        vals = vals(:);
        vals = vals(isfinite(vals));
        x = plot_r + offsets(c);

        if isempty(vals)
            continue
        end

        boxchart(ax, x * ones(size(vals)), vals, ...
            'BoxFaceColor', colors(c,:), ...
            'BoxFaceAlpha', 0.55, ...
            'MarkerStyle', 'none', ...
            'BoxWidth', 0.26, ...
            'LineWidth', 0.45);

        if numel(vals) > max_scatter
            vals_for_scatter = vals(randperm(numel(vals), max_scatter));
        else
            vals_for_scatter = vals;
        end

        scatter(ax, x + 0.035*randn(size(vals_for_scatter)), vals_for_scatter, ...
            8, colors(c,:), 'filled', 'MarkerFaceAlpha', 0.14, ...
            'MarkerEdgeColor', 'none');

        all_plot_vals = [all_plot_vals; vals]; %#ok<AGROW>
    end
end

if isempty(all_plot_vals)
    warning('No %s values available for combined compact summary.', feature_name);
    hold(ax, 'off');
    return
end

[y_min, y_max, y_range] = set_zoomed_ylim(all_plot_vals, 2.7);
y_limits = [y_min, y_max];
y_pad = 0.08 * max(y_range, eps);

for plot_r = 1:n_regions_to_plot
    r = region_inds_to_plot(plot_r);
    y_text = y_limits(2) + 0.38*y_pad;
    stats_txt = sprintf('Welch''s t(%d)=%.2f\np=%.2g, d=%.2f', ...
        stats_allpairs.tDf(r), ...
        stats_allpairs.tStat(r), ...
        stats_allpairs.pValue(r), ...
        stats_allpairs.CohensD(r));
    text(ax, plot_r, y_text, stats_txt, 'HorizontalAlignment', 'center', 'FontSize', 8);

    if is_significant_allpairs(stats_allpairs, r)
        x1 = plot_r + offsets(1);
        x2 = plot_r + offsets(2);
        y_bracket = y_limits(2) + 1.35*y_pad;
        y_tick = 0.12*y_pad;
        plot(ax, [x1 x1 x2 x2], ...
            [y_bracket-y_tick, y_bracket, y_bracket, y_bracket-y_tick], ...
            'k-', 'LineWidth', 0.5);
        text(ax, plot_r, y_bracket + 0.08*y_pad, p_to_star(stats_allpairs.pValue(r)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 13, ...
            'FontWeight', 'bold');
    end
end

xlim(ax, [0.45, n_regions_to_plot + 0.55]);
xticks(ax, 1:n_regions_to_plot);
xticklabels(ax, region_labels_to_plot);
xtickangle(ax, 25);
ylabel(ax, [feature_name, ' (bits)']);
title(ax, feature_name, 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'off');
set(ax, 'FontSize', 11, 'LineWidth', 0.45);

if show_legend
    legend_handles = gobjects(2,1);
    for c = 1:2
        legend_handles(c) = scatter(ax, nan, nan, 30, colors(c,:), 'filled', 'MarkerEdgeColor', 'none');
    end
    legend(ax, legend_handles, activity_labels(1:2), 'Location', 'northoutside', ...
        'Orientation', 'horizontal', 'Box', 'off');
end
hold(ax, 'off');
end

function print_editable_svg(fig_handle, svg_file)
% Use MATLAB's vector renderer so Illustrator receives editable paths/text
% instead of a raster snapshot embedded inside an SVG wrapper.
set(fig_handle, 'Renderer', 'painters');
print(fig_handle, svg_file, '-dsvg', '-vector');
end

function [y_min, y_max, y_range] = set_zoomed_ylim(all_plot_vals, top_padding_units)
% Robust y-limits for paper summary plots: preserve the statistics/data, but
% avoid rare extreme CoI values flattening the visible distributions.
finite_vals = all_plot_vals(isfinite(all_plot_vals));
if isempty(finite_vals)
    y_min = 0;
    y_max = 1;
    y_range = 1;
    ylim([0 1]);
    return
end

if numel(finite_vals) >= 20
    y_min = prctile(finite_vals, 1);
    y_max = prctile(finite_vals, 99);
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

y_pad = 0.08 * max(y_range, eps);
ylim([y_min - y_pad, y_max + top_padding_units*y_pad]);
end

function is_sig = is_significant_allpairs(stats_allpairs, row_idx)
% Prefer the FDR decision when it has already been added to the all-pairs
% table; otherwise fall back to the raw Welch p-value.
if ismember('FDR_sig', stats_allpairs.Properties.VariableNames)
    is_sig = logical(stats_allpairs.FDR_sig(row_idx));
else
    is_sig = stats_allpairs.pValue(row_idx) < 0.05;
end
end

function star_txt = p_to_star(p_val)
if ~isfinite(p_val)
    star_txt = '';
elseif p_val < 0.001
    star_txt = '***';
elseif p_val < 0.01
    star_txt = '**';
else
    star_txt = '*';
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
comparison_tag = matlab.lang.makeValidName(strjoin(activity_tags(1:2), '_'));
end

function difference_tag = make_activity_difference_tag(activity_tags)
difference_tag = matlab.lang.makeValidName([activity_tags{2} '_minus_' activity_tags{1}]);
end

function CoIData = get_cached_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, times)
% Avoid reloading the same participant/activity CoI files for every region
% and for both redundancy/synergy passes. The raw CoIData object contains
% all regions and both masks, so it is safe to reuse.

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
% Simple Benjamini-Hochberg FDR correction returning adjusted p-values and significance mask
% pvals: vector of p-values
% alpha: target FDR level (e.g., 0.05)
adj_p = nan(size(pvals));
sig = false(size(pvals));
valid = ~isnan(pvals);
if ~any(valid)
    return;
end
pv = pvals(valid);
[sorted_p, sort_idx] = sort(pv);
rank = (1:numel(sorted_p))';
% BH adjusted p (store min for increasing rank)
adj_sorted = sorted_p .* numel(sorted_p) ./ rank;
% enforce monotonicity
for i = numel(adj_sorted)-1:-1:1
    adj_sorted(i) = min(adj_sorted(i), adj_sorted(i+1));
end
% cap at 1
adj_sorted(adj_sorted>1) = 1;
% place back
adj_p(valid) = adj_sorted(sort_idx);
% determine significance
sig(valid) = adj_p(valid) <= alpha;
end
