function compare_region_patterns(participants, basefold, datatype, activity_tag, condition, cutoff, timing, out_dir)
%COMPARE_REGION_PATTERNS Compare Co-I / redundancy / synergy patterns across regions
%
% This function performs a group-level comparison of time-by-time
% co-information (Co-I) matrices across anatomical region groups.
% It is intended to be run after each fly's pairwise Co-I matrices
% have been computed and packaged by `max_get_plotting_CoI_raw_rf`.
%
% High-level steps:
% 1. Aggregate per-pair Co-I matrices across flies for each region.
% 2. Compute group mean Co-I matrix for each region and split into
%    redundancy (positive values) and synergy (negative values).
% 3. For every pair of regions, compute similarity metrics between
%    their mean matrices: Pearson correlation (r) and SSIM (structural
%    similarity index). Similarity is computed separately for
%    redundancy and synergy components.
% 4. Assess significance with two permutation nulls:
%    - circular time-shift (preserves temporal autocorrelation)
%    - entry shuffle (global scramble of values)
%    For both Pearson r and SSIM we compute null distributions
%    and p-values (one-sided: observed >= null).
% 5. Compute bootstrap CIs for metrics by resampling pairwise
%    matrices (across pairs/fly contributions).
% 6. Save results as a CSV and generate summary figures (heatmaps,
%    null histograms, observed - null differences, z-scores).
%
% Notes on conventions used here:
% - Redundancy = positive part of the mean Co-I matrix
% - Synergy = absolute value of the negative part (so shapes are
%   compared on a positive scale)
% - Masking: comparisons are restricted to a union of significant
%   cells from the two regions by default; fallback is all finite cells.
% - SSIM requires the Image Processing Toolbox; if installed, true
%   SSIM is used (we tested with exist('ssim','file') earlier).
%
% Usage:
% compare_region_patterns(participants, basefold, datatype, activity_tag, condition, cutoff, timing, out_dir)
%
% Inputs:
% - participants: cell array of subject names (same as used elsewhere)
% - basefold, datatype, activity_tag, condition, cutoff: passed to
%   max_get_plotting_CoI_raw_rf (same as in plot_all_flies_average_2.m)
% - timing: vector of time points (e.g., -25:74)
% - out_dir: directory where outputs (figures, CSV) will be saved
%
% The function computes group mean matrices per region, splits
% redundant (positive) and synergetic (negative) components, then
% computes pairwise similarity (Pearson r) between regions for
% redundancy and synergy separately. It runs permutation tests (circular
% time-shift and entry-shuffle) and bootstraps CIs across pairs.

if nargin < 9 || isempty(out_dir)
    out_dir = fullfile(basefold, 'Plots', 'Region_Patterns');
end
if ~exist(out_dir,'dir'), mkdir(out_dir); end

wb = [];
try
    wb = waitbar(0, 'Starting region comparison analysis...');
catch
    wb = [];
end
cleanup_wb = onCleanup(@()close_waitbar_safe(wb));

regions = {'peripheral','central','central_peripheral','peripheral_peripheral','central_central'};
region_labels = {'Peripheral','Central','Central-Peripheral','Peripheral-Peripheral','Central-Central'};
n_regions = numel(regions);

% -----------------------------------------------------------------
% Aggregate across flies: collect all pairwise Co-I arrays and masks
% for each anatomical region. The resulting `group` struct contains
% 3D arrays (time x time x n_pairs_total) and mask arrays so we can
% compute group-level means and estimate variability.
% -----------------------------------------------------------------
% Aggregate across flies (reuse same logic as plot_all_flies_average_2)
for r = 1:n_regions
    group.(regions{r}).FFi_all = [];
    group.(regions{r}).mask_all = [];
    group.(regions{r}).mask_r_all = [];
    group.(regions{r}).mask_s_all = [];
    group.(regions{r}).mi1_all = [];
    group.(regions{r}).mi2_all = [];
    group.(regions{r}).n_pairs = 0;
end

% Loop over participants and append each fly's pairwise results into
% the region-specific accumulators. Each `fly.FFi_all` is expected
% to be time x time x n_pairs_for_this_fly. We concatenate along the
% 3rd dim to produce a pooled set of pair matrices.
for s = 1:numel(participants)
    participant_name = char(participants(s));
    % Load this fly's raw CoI struct (fields per region)
    flyData = max_get_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, timing);
    for r = 1:n_regions
        reg = regions{r};
        fly = flyData.(reg);
        if isempty(fly.FFi_all)
            % No pairs for this region in this fly
            continue;
        end
        % Concatenate: time x time x (pairs accumulated so far)
        group.(reg).FFi_all = cat(3, group.(reg).FFi_all, fly.FFi_all);
        % Collect masks and MI vectors for downstream bootstrapping
        group.(reg).mask_all   = cat(3, group.(reg).mask_all, fly.mask);
        group.(reg).mask_r_all = cat(3, group.(reg).mask_r_all, fly.mask_r);
        group.(reg).mask_s_all = cat(3, group.(reg).mask_s_all, fly.mask_s);
        group.(reg).mi1_all = [group.(reg).mi1_all, fly.mi1];
        group.(reg).mi2_all = [group.(reg).mi2_all, fly.mi2];
        group.(reg).n_pairs = group.(reg).n_pairs + size(fly.FFi_all,3);
    end
    if ~isempty(wb)
        waitbar(0.10 * s / max(1, numel(participants)), wb, sprintf('Aggregating fly %d/%d...', s, numel(participants)));
    end
end

% -----------------------------------------------------------------
% Compute group mean Co-I per region and split into redundancy and
% synergy components. We also compute an aggregate significance mask
% (proportion > 0 across pairwise masks) to focus comparisons on
% cells that are reliably non-null across pairs/flies.
% -----------------------------------------------------------------
% Compute mean and split redundant / synergetic, keep masks
for r = 1:n_regions
    reg = regions{r};
    if isempty(group.(reg).FFi_all)
        warning('Region %s empty, skipping', reg);
        region_mean = nan(length(timing));
        region_mask = false(length(timing));
        region_n = 0;
    else
        region_mean = nanmean(group.(reg).FFi_all, 3);
        region_mask = mean(group.(reg).mask_all, 3, 'omitnan') > 0; % proportion >0
        region_n = group.(reg).n_pairs;
    end
    region_struct(r).name = reg;
    region_struct(r).label = region_labels{r};
    region_struct(r).mean = region_mean;
    region_struct(r).redundant = region_mean .* (region_mean > 0);
    region_struct(r).synergetic = region_mean .* (region_mean < 0);
    region_struct(r).mask = region_mask;
    region_struct(r).n = region_n;
    % keep the raw stacked pairwise arrays too (used by bootstrapping)
    region_struct(r).FFi_all = group.(reg).FFi_all;
end

% -----------------------------------------------------------------
% Prepare outputs and storage for pairwise statistics
% -----------------------------------------------------------------
% Prepare outputs
Nperm = 1000; Nboot = 1000; % defaults
metrics = {'Pearson','SSIM'};

% Preallocate result tables
pairs = nchoosek(1:n_regions,2);
n_pairs = size(pairs,1);
results = table;
results.RegionA = strings(n_pairs,1);
results.RegionB = strings(n_pairs,1);
results.r_red = nan(n_pairs,1);
results.p_circ_red = nan(n_pairs,1);
results.p_entry_red = nan(n_pairs,1);
results.CI_low_red = nan(n_pairs,1);
results.CI_high_red = nan(n_pairs,1);
results.ssim_red = nan(n_pairs,1);

results.r_syn = nan(n_pairs,1);
results.p_circ_syn = nan(n_pairs,1);
results.p_entry_syn = nan(n_pairs,1);
results.CI_low_syn = nan(n_pairs,1);
results.CI_high_syn = nan(n_pairs,1);
results.ssim_syn = nan(n_pairs,1);

% SSIM permutation results columns (will be written to CSV)
results.ssim_p_circ_red = nan(n_pairs,1);
results.ssim_p_entry_red = nan(n_pairs,1);
results.ssim_CI_low_red = nan(n_pairs,1);
results.ssim_CI_high_red = nan(n_pairs,1);

results.ssim_p_circ_syn = nan(n_pairs,1);
results.ssim_p_entry_syn = nan(n_pairs,1);
results.ssim_CI_low_syn = nan(n_pairs,1);
results.ssim_CI_high_syn = nan(n_pairs,1);

% Signed (whole mean matrix) SSIM: observed, permutation p-values and CIs
results.ssim_signed = nan(n_pairs,1);
results.ssim_p_circ_signed = nan(n_pairs,1);
results.ssim_p_entry_signed = nan(n_pairs,1);
results.ssim_CI_low_signed = nan(n_pairs,1);
results.ssim_CI_high_signed = nan(n_pairs,1);

% Prepare storage for null distributions
nulls = struct('circ_red', cell(n_pairs,1), 'entry_red', cell(n_pairs,1), 'circ_syn', cell(n_pairs,1), 'entry_syn', cell(n_pairs,1), ...
    'circ_ssim_red', cell(n_pairs,1), 'entry_ssim_red', cell(n_pairs,1), 'circ_ssim_syn', cell(n_pairs,1), 'entry_ssim_syn', cell(n_pairs,1), ...
    'circ_ssim_signed', cell(n_pairs,1), 'entry_ssim_signed', cell(n_pairs,1));

% -----------------------------------------------------------------
% Loop over all unordered region pairs and compute similarity metrics.
% For each pair we:
% - select a mask (union of significant cells)
% - compute Pearson r and SSIM for redundancy and synergy
% - compute permutation p-values (circular shift and entry shuffle)
% - compute bootstrap CIs
% - store observed values, p-values and null distributions
% -----------------------------------------------------------------
% Loop over pairs
for k = 1:n_pairs
    i = pairs(k,1); j = pairs(k,2);
    A = region_struct(i); B = region_struct(j);
    results.RegionA(k) = A.label;
    results.RegionB(k) = B.label;

    % Use union mask by default: include cells significant in either
    % region so comparisons are sensitive to shared structure.
    mask_union = (A.mask > 0) | (B.mask > 0);
    if all(~mask_union(:))
        % fallback to all non-NaN
        mask_union = ~isnan(A.mean) & ~isnan(B.mean);
    end

    % ------------------- Redundancy comparisons -------------------
    % Use the positive part of the mean matrix (redundancy).
    % Keep zeros for non-redundant cells so spatial pattern is retained.
    Xr = A.redundant; Yr = B.redundant;
    % Pearson r nulls
    [r_obs_red, p_circ_red, r_null_red] = perm_test_circshift(Xr, Yr, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.02) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: redundancy Pearson', k, n_pairs)); end
    [~, p_entry_red, r_null_entry_red] = perm_test_entry_shuffle(Xr, Yr, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.05) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: redundancy shuffle', k, n_pairs)); end
    [r_boot_mean_red, CI_red] = bootstrap_similarity(A.FFi_all .* (A.mean>0), B.FFi_all .* (B.mean>0), mask_union, Nboot);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.08) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: redundancy bootstrap', k, n_pairs)); end
    % SSIM observed and nulls
    ssim_obs_red = matrix_ssim(Xr, Yr, mask_union);
    [ssim_p_circ_red, ~, ssim_null_circ_red] = perm_test_circshift_metric(@matrix_ssim, Xr, Yr, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.11) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: redundancy SSIM', k, n_pairs)); end
    [~, ssim_p_entry_red, ssim_null_entry_red] = perm_test_entry_shuffle_metric(@matrix_ssim, Xr, Yr, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.14) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: redundancy SSIM shuffle', k, n_pairs)); end
    % bootstrap CI for SSIM
    [ssim_boot_mean_red, ssim_CI_red] = bootstrap_ssim(A.FFi_all .* (A.mean>0), B.FFi_all .* (B.mean>0), mask_union, Nboot);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.17) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: redundancy SSIM bootstrap', k, n_pairs)); end
    results.r_red(k) = r_obs_red;
    results.p_circ_red(k) = p_circ_red;
    results.p_entry_red(k) = p_entry_red;
    results.CI_low_red(k) = CI_red(1);
    results.CI_high_red(k) = CI_red(2);
    results.ssim_p_circ_red(k) = ssim_p_circ_red;
    results.ssim_p_entry_red(k) = ssim_p_entry_red;
    results.ssim_CI_low_red(k) = ssim_CI_red(1);
    results.ssim_CI_high_red(k) = ssim_CI_red(2);
    % Store SSIM observed value (matrix_ssim normalizes and calls
    % MATLAB's ssim()). Use try/catch to be robust to edge cases.
    try
        results.ssim_red(k) = matrix_ssim(Xr, Yr, mask_union);
    catch
        results.ssim_red(k) = NaN;
    end

    % ------------------- Synergy comparisons ---------------------
    % For synergy we compare the absolute value of negative entries so
    % that synergetic 'hotspots' (negative values) are compared as
    % positive structures.
    Xs = abs(A.synergetic); Ys = abs(B.synergetic);
    [r_obs_syn, p_circ_syn, r_null_syn] = perm_test_circshift(Xs, Ys, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.20) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: synergy Pearson', k, n_pairs)); end
    [~, p_entry_syn, r_null_entry_syn] = perm_test_entry_shuffle(Xs, Ys, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.23) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: synergy shuffle', k, n_pairs)); end
    [r_boot_mean_syn, CI_syn] = bootstrap_similarity(abs(A.FFi_all .* (A.mean<0)), abs(B.FFi_all .* (B.mean<0)), mask_union, Nboot);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.26) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: synergy bootstrap', k, n_pairs)); end
    % SSIM observed and nulls for synergy
    ssim_obs_syn = matrix_ssim(Xs, Ys, mask_union);
    [ssim_p_circ_syn, ~, ssim_null_circ_syn] = perm_test_circshift_metric(@matrix_ssim, Xs, Ys, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.29) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: synergy SSIM', k, n_pairs)); end
    [~, ssim_p_entry_syn, ssim_null_entry_syn] = perm_test_entry_shuffle_metric(@matrix_ssim, Xs, Ys, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.32) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: synergy SSIM shuffle', k, n_pairs)); end
    [ssim_boot_mean_syn, ssim_CI_syn] = bootstrap_ssim(abs(A.FFi_all .* (A.mean<0)), abs(B.FFi_all .* (B.mean<0)), mask_union, Nboot);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.35) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: synergy SSIM bootstrap', k, n_pairs)); end
    results.r_syn(k) = r_obs_syn;
    results.p_circ_syn(k) = p_circ_syn;
    results.p_entry_syn(k) = p_entry_syn;
    results.CI_low_syn(k) = CI_syn(1);
    results.CI_high_syn(k) = CI_syn(2);
    results.ssim_p_circ_syn(k) = ssim_p_circ_syn;
    results.ssim_p_entry_syn(k) = ssim_p_entry_syn;
    results.ssim_CI_low_syn(k) = ssim_CI_syn(1);
    results.ssim_CI_high_syn(k) = ssim_CI_syn(2);
    try
        results.ssim_syn(k) = matrix_ssim(Xs, Ys, mask_union);
    catch
        results.ssim_syn(k) = NaN;
    end
    % ------------------- Signed whole-matrix SSIM -------------------
    % Compute SSIM on the signed group mean matrices (no sign-splitting)
    Xsld = A.mean; Ysld = B.mean;
    try
        ssim_obs_signed = matrix_ssim(Xsld, Ysld, mask_union);
    catch
        ssim_obs_signed = NaN;
    end
    [ssim_p_circ_signed, ~, ssim_null_circ_signed] = perm_test_circshift_metric(@matrix_ssim, Xsld, Ysld, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.38) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: signed SSIM', k, n_pairs)); end
    [~, ssim_p_entry_signed, ssim_null_entry_signed] = perm_test_entry_shuffle_metric(@matrix_ssim, Xsld, Ysld, mask_union, Nperm);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.41) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: signed SSIM shuffle', k, n_pairs)); end
    [ssim_boot_mean_signed, ssim_CI_signed] = bootstrap_ssim(A.FFi_all, B.FFi_all, mask_union, Nboot);
    if ~isempty(wb), waitbar(0.10 + ((k-1) + 0.44) / max(1,n_pairs) * 0.80, wb, sprintf('Pair %d/%d: signed SSIM bootstrap', k, n_pairs)); end
    results.ssim_signed(k) = ssim_obs_signed;
    results.ssim_p_circ_signed(k) = ssim_p_circ_signed;
    results.ssim_p_entry_signed(k) = ssim_p_entry_signed;
    results.ssim_CI_low_signed(k) = ssim_CI_signed(1);
    results.ssim_CI_high_signed(k) = ssim_CI_signed(2);
    % store null distributions for plotting later (r and SSIM)
    nulls(k).circ_red = r_null_red;
    nulls(k).entry_red = r_null_entry_red;
    nulls(k).circ_syn = r_null_syn;
    nulls(k).entry_syn = r_null_entry_syn;
    % add SSIM nulls
    nulls(k).circ_ssim_red = ssim_null_circ_red;
    nulls(k).entry_ssim_red = ssim_null_entry_red;
    nulls(k).circ_ssim_syn = ssim_null_circ_syn;
    nulls(k).entry_ssim_syn = ssim_null_entry_syn;
    % signed SSIM nulls
    nulls(k).circ_ssim_signed = ssim_null_circ_signed;
    nulls(k).entry_ssim_signed = ssim_null_entry_signed;
end

% -----------------------------------------------------------------
% Apply Benjamini-Hochberg FDR correction to a selected subset of p-values
% We only correct: circular-shift p-values for Pearson r (redundancy/synergy)
% and circular-shift SSIM p-values for the signed mean matrices.
% -----------------------------------------------------------------
% Columns to correct
pval_cols = {'p_circ_red', 'p_circ_syn', 'ssim_p_circ_signed'};

% Collect p-values and masks for mapping back
all_pvals = [];
col_sizes = zeros(numel(pval_cols),1);
for ci = 1:numel(pval_cols)
    p = results.(pval_cols{ci});
    valid_mask = ~isnan(p);
    col_sizes(ci) = sum(valid_mask);
    all_pvals = [all_pvals; p(valid_mask)];
end

% Benjamini-Hochberg adjustment (standard step-up procedure producing adjusted p-values)
if ~isempty(all_pvals)
    m = numel(all_pvals);
    [p_sorted, sort_idx] = sort(all_pvals, 'ascend');
    adj_sorted = nan(m,1);
    % Compute p_adj_sorted(i) = min_{j>=i} (p_sorted(j) * m / j)
    for i = m:-1:1
        adj_val = p_sorted(i) * m / i;
        if i < m
            adj_sorted(i) = min(adj_val, adj_sorted(i+1));
        else
            adj_sorted(i) = adj_val;
        end
    end
    adj_sorted = min(adj_sorted, 1);
    % map back to original order
    adj_all = nan(m,1);
    adj_all(sort_idx) = adj_sorted;

    % distribute corrected p-values back into results table columns
    ptr = 1;
    for ci = 1:numel(pval_cols)
        col_name = pval_cols{ci};
        p = results.(col_name);
        valid_mask = ~isnan(p);
        n = col_sizes(ci);
        if n > 0
            corrected_vals = adj_all(ptr:ptr+n-1);
            p(valid_mask) = corrected_vals;
            results.(col_name) = p;
            ptr = ptr + n;
        end
    end
end

% Save numeric summary table to CSV for later inspection and reporting.
% NOTE: Circular-shift p-values for Pearson r (redundancy & synergy) and
% signed-SSIM have been corrected using Benjamini-Hochberg FDR (alpha = 0.05).
outcsv = fullfile(out_dir, sprintf('region_pattern_similarity_%s_%s.csv', activity_tag, condition));
writetable(results, outcsv);

% -----------------------------------------------------------------
% Make summary heatmaps: Pearson r for redundancy and synergy
% NOTE: Significance stars (*) indicate p < 0.05 after Benjamini-Hochberg
% FDR correction applied to circular-shift p-values for Pearson r and
% signed-SSIM (only these p-values were corrected).
% -----------------------------------------------------------------
% Make simple heatmaps for r_red and r_syn (regions x regions)
R_red = nan(n_regions); R_syn = nan(n_regions);
P_red = ones(n_regions); P_syn = ones(n_regions);
for k = 1:n_pairs
    i = pairs(k,1); j = pairs(k,2);
    R_red(i,j) = results.r_red(k); R_red(j,i) = R_red(i,j);
    R_syn(i,j) = results.r_syn(k); R_syn(j,i) = R_syn(i,j);
    P_red(i,j) = results.p_circ_red(k); P_red(j,i) = P_red(i,j);
    P_syn(i,j) = results.p_circ_syn(k); P_syn(j,i) = P_syn(i,j);
end
for d = 1:n_regions
    R_red(d,d) = 1; R_syn(d,d) = 1; P_red(d,d) = 0; P_syn(d,d)=0;
end

fig = figure('Position',[100 100 1000 420]);
subplot(1,2,1);
imagesc(R_red, [-1 1]); colormap(gca, parula); colorbar; axis square;
title('Redundant pattern similarity (r)');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);
% annotate significance
hold on;
sig = P_red < 0.05;
[ri,cj] = find(sig);
for t = 1:numel(ri)
    text(cj(t), ri(t), '*','HorizontalAlignment','center','Color','k','FontSize',14);
end
hold off;

subplot(1,2,2);
imagesc(R_syn, [-1 1]); colormap(gca, parula); colorbar; axis square;
title('Synergetic pattern similarity (r)');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);
hold on;
sig2 = P_syn < 0.05;
[ri,cj] = find(sig2);
for t = 1:numel(ri)
    text(cj(t), ri(t), '*','HorizontalAlignment','center','Color','k','FontSize',14);
end
hold off;

% Persist the main Pearson heatmaps
saveas(fig, fullfile(out_dir, sprintf('region_similarity_%s_%s.png', activity_tag, condition)));
print(fig, fullfile(out_dir, sprintf('region_similarity_%s_%s.svg', activity_tag, condition)), '-dsvg');
close(fig);

% -----------------------------------------------------------------
% SSIM heatmaps: complementary view that emphasizes local structural
% similarity rather than linear cell-by-cell correspondence.
% -----------------------------------------------------------------
% --- SSIM heatmaps ---
S_red = nan(n_regions); S_syn = nan(n_regions);
for k = 1:n_pairs
    i = pairs(k,1); j = pairs(k,2);
    S_red(i,j) = results.ssim_red(k); S_red(j,i) = S_red(i,j);
    S_syn(i,j) = results.ssim_syn(k); S_syn(j,i) = S_syn(i,j);
end
for d = 1:n_regions
    S_red(d,d) = 1; S_syn(d,d) = 1;
end
figSmap = figure('Position',[100 100 1000 420]);
subplot(1,2,1);
imagesc(S_red, [0 1]); colormap(gca, parula); colorbar; axis square;
title('Redundant pattern similarity (SSIM)');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);

subplot(1,2,2);
imagesc(S_syn, [0 1]); colormap(gca, parula); colorbar; axis square;
title('Synergetic pattern similarity (SSIM)');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);
% Annotate SSIM heatmaps with permutation significance (entry-shuffle stars
% and circular-shift p-values). Build p-value matrices from results.
P_ssim_entry_red = ones(n_regions);
P_ssim_circ_red  = ones(n_regions);
P_ssim_entry_syn = ones(n_regions);
P_ssim_circ_syn  = ones(n_regions);
for k = 1:n_pairs
    i = pairs(k,1); j = pairs(k,2);
    P_ssim_entry_red(i,j) = results.ssim_p_entry_red(k); P_ssim_entry_red(j,i) = P_ssim_entry_red(i,j);
    P_ssim_circ_red(i,j)  = results.ssim_p_circ_red(k);  P_ssim_circ_red(j,i)  = P_ssim_circ_red(i,j);
    P_ssim_entry_syn(i,j) = results.ssim_p_entry_syn(k); P_ssim_entry_syn(j,i) = P_ssim_entry_syn(i,j);
    P_ssim_circ_syn(i,j)  = results.ssim_p_circ_syn(k);  P_ssim_circ_syn(j,i)  = P_ssim_circ_syn(i,j);
end
for d = 1:n_regions
    P_ssim_entry_red(d,d) = NaN; P_ssim_circ_red(d,d) = NaN;
    P_ssim_entry_syn(d,d) = NaN; P_ssim_circ_syn(d,d) = NaN;
end

% Now overlay annotations on the SSIM heatmaps
% We'll use circular-shift p-values for star significance and overlay the
% exact SSIM numeric value in each cell for clarity.
% Left: redundancy
subplot(1,2,1);
hold on;
for rrow = 1:n_regions
    for ccol = 1:n_regions
        sval = S_red(rrow, ccol);
        if ~isnan(sval)
            % SSIM value (centered)
            text(ccol, rrow, sprintf('%.3f', sval), 'HorizontalAlignment','center', 'Color','k', 'FontSize',10);
        end
        p_circ = P_ssim_circ_red(rrow, ccol);
        % star for circular-shift significance
        if ~isnan(p_circ) && p_circ < 0.05
            text(ccol, rrow-0.25, '*', 'HorizontalAlignment','center', 'Color','k', 'FontSize',14);
        end
    end
end
hold off;

% Right: synergy
subplot(1,2,2);
hold on;
for rrow = 1:n_regions
    for ccol = 1:n_regions
        sval = S_syn(rrow, ccol);
        if ~isnan(sval)
            text(ccol, rrow, sprintf('%.3f', sval), 'HorizontalAlignment','center', 'Color','k', 'FontSize',10);
        end
        p_circ = P_ssim_circ_syn(rrow, ccol);
        if ~isnan(p_circ) && p_circ < 0.05
            text(ccol, rrow-0.25, '*', 'HorizontalAlignment','center', 'Color','k', 'FontSize',14);
        end
    end
end
hold off;

% Save annotated SSIM figures
saveas(figSmap, fullfile(out_dir, sprintf('region_similarity_ssim_annotated_%s_%s.png', activity_tag, condition)));
print(figSmap, fullfile(out_dir, sprintf('region_similarity_ssim_annotated_%s_%s.svg', activity_tag, condition)), '-dsvg');
close(figSmap);

% --- Signed SSIM heatmap (whole signed mean matrices) ---
S_signed = nan(n_regions);
P_ssim_circ_signed = ones(n_regions);
for k = 1:n_pairs
    i = pairs(k,1); j = pairs(k,2);
    S_signed(i,j) = results.ssim_signed(k); S_signed(j,i) = S_signed(i,j);
    P_ssim_circ_signed(i,j) = results.ssim_p_circ_signed(k); P_ssim_circ_signed(j,i) = P_ssim_circ_signed(i,j);
end
for d = 1:n_regions
    if isnan(S_signed(d,d)), S_signed(d,d) = 1; end
    P_ssim_circ_signed(d,d) = NaN;
end

figSsig = figure('Position',[100 100 600 500]);
imagesc(S_signed);
colormap(parula); colorbar; axis square;
title('Signed mean Co-I similarity (SSIM)');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);
% overlay numeric SSIM and circular-shift significance stars
hold on;
for rrow = 1:n_regions
    for ccol = 1:n_regions
        sval = S_signed(rrow, ccol);
        if ~isnan(sval)
            text(ccol, rrow, sprintf('%.3f', sval), 'HorizontalAlignment','center', 'Color','k', 'FontSize',10);
        end
        p_circ = P_ssim_circ_signed(rrow, ccol);
        if ~isnan(p_circ) && p_circ < 0.05
            text(ccol, rrow-0.25, '*', 'HorizontalAlignment','center', 'Color','k', 'FontSize',14);
        end
    end
end
hold off;

saveas(figSsig, fullfile(out_dir, sprintf('region_similarity_ssim_signed_annotated_%s_%s.png', activity_tag, condition)));
print(figSsig, fullfile(out_dir, sprintf('region_similarity_ssim_signed_annotated_%s_%s.svg', activity_tag, condition)), '-dsvg');
close(figSsig);

if ~isempty(wb)
    waitbar(0.95, wb, 'Computing final summaries and writing outputs...');
end

% --- Condition vs null: compute mean and std of circular-shift nulls and plot difference/z-score ---
MeanNull_red = nan(n_regions); StdNull_red = nan(n_regions);
MeanNull_syn = nan(n_regions); StdNull_syn = nan(n_regions);
for k = 1:n_pairs
    i = pairs(k,1); j = pairs(k,2);
    nr = nulls(k).circ_red;
    ns = nulls(k).circ_syn;
    MeanNull_red(i,j) = mean(nr); MeanNull_red(j,i) = MeanNull_red(i,j);
    StdNull_red(i,j) = std(nr); StdNull_red(j,i) = StdNull_red(i,j);
    MeanNull_syn(i,j) = mean(ns); MeanNull_syn(j,i) = MeanNull_syn(i,j);
    StdNull_syn(i,j) = std(ns); StdNull_syn(j,i) = StdNull_syn(i,j);
end
for d = 1:n_regions
    MeanNull_red(d,d) = 0; StdNull_red(d,d) = NaN;
    MeanNull_syn(d,d) = 0; StdNull_syn(d,d) = NaN;
end

Diff_red = R_red - MeanNull_red; % observed minus average null
Z_red = (R_red - MeanNull_red) ./ StdNull_red; % z-score relative to circular null
Diff_syn = R_syn - MeanNull_syn;
Z_syn = (R_syn - MeanNull_syn) ./ StdNull_syn;

% Plot difference heatmaps (observed minus mean null)
figD = figure('Position',[100 100 1000 420]);
subplot(1,2,1);
imagesc(Diff_red, [-0.5 0.5]); colormap(gca, parula); colorbar; axis square;
title('Redundancy: observed r - mean(null_{circ})');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);

subplot(1,2,2);
imagesc(Diff_syn, [-0.5 0.5]); colormap(gca, parula); colorbar; axis square;
title('Synergy: observed r - mean(null_{circ})');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);

saveas(figD, fullfile(out_dir, sprintf('region_similarity_diff_null_%s_%s.png', activity_tag, condition)));
print(figD, fullfile(out_dir, sprintf('region_similarity_diff_null_%s_%s.svg', activity_tag, condition)), '-dsvg');
close(figD);

% Plot z-score heatmaps
figZ = figure('Position',[100 100 1000 420]);
subplot(1,2,1);
imagesc(Z_red, [-3 3]); colormap(gca, parula); colorbar; axis square;
title('Redundancy: z-score vs circular null');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);

subplot(1,2,2);
imagesc(Z_syn, [-3 3]); colormap(gca, parula); colorbar; axis square;
title('Synergy: z-score vs circular null');
xticks(1:n_regions); yticks(1:n_regions); xticklabels(region_labels); yticklabels(region_labels);

saveas(figZ, fullfile(out_dir, sprintf('region_similarity_zscore_null_%s_%s.png', activity_tag, condition)));
print(figZ, fullfile(out_dir, sprintf('region_similarity_zscore_null_%s_%s.svg', activity_tag, condition)), '-dsvg');
close(figZ);

% --- Comparison figures: observed vs permutation nulls ---
% Redundancy null comparison
figR = figure('Name','Redundancy nulls','Units','pixels','Position',[100 100 1200 800]);
rows = ceil(n_pairs/2); cols = 4; % two columns per pair (circ, entry) displayed side-by-side
tl = tiledlayout(rows, cols, 'TileSpacing','Compact','Padding','Compact');
title(tl, sprintf('Redundancy: observed r vs null distributions (%s, %s)', activity_tag, condition));
for k = 1:n_pairs
    % determine tile indexes to place circshift and entry histograms
    tile_idx = (k-1)*2 + 1;
    nexttile(tile_idx);
    h1 = histogram(nulls(k).circ_red, 'Normalization','probability','FaceColor',[0.8 0.8 0.8]); hold on;
    yl = ylim();
    plot([results.r_red(k) results.r_red(k)], [0 yl(2)], 'k-', 'LineWidth',2);
    title(sprintf('%s vs %s (circ)', results.RegionA(k), results.RegionB(k)));
    xlabel('r'); ylabel('prob'); hold off;

    nexttile(tile_idx+1);
    histogram(nulls(k).entry_red, 'Normalization','probability','FaceColor',[0.7 0.7 0.7]); hold on;
    yl = ylim();
    plot([results.r_red(k) results.r_red(k)], [0 yl(2)], 'k-', 'LineWidth',2);
    title(sprintf('%s vs %s (entry)', results.RegionA(k), results.RegionB(k)));
    xlabel('r'); ylabel('prob'); hold off;
end
saveas(figR, fullfile(out_dir, sprintf('redundancy_nulls_%s_%s.png', activity_tag, condition)));
print(figR, fullfile(out_dir, sprintf('redundancy_nulls_%s_%s.svg', activity_tag, condition)), '-dsvg');
close(figR);

% Synergy null comparison
figS = figure('Name','Synergy nulls','Units','pixels','Position',[100 100 1200 800]);
tl2 = tiledlayout(rows, cols, 'TileSpacing','Compact','Padding','Compact');
title(tl2, sprintf('Synergy: observed r vs null distributions (%s, %s)', activity_tag, condition));
for k = 1:n_pairs
    tile_idx = (k-1)*2 + 1;
    nexttile(tile_idx);
    histogram(nulls(k).circ_syn, 'Normalization','probability','FaceColor',[0.8 0.8 0.8]); hold on;
    yl = ylim();
    plot([results.r_syn(k) results.r_syn(k)], [0 yl(2)], 'k-', 'LineWidth',2);
    title(sprintf('%s vs %s (circ)', results.RegionA(k), results.RegionB(k)));
    xlabel('r'); ylabel('prob'); hold off;

    nexttile(tile_idx+1);
    histogram(nulls(k).entry_syn, 'Normalization','probability','FaceColor',[0.7 0.7 0.7]); hold on;
    yl = ylim();
    plot([results.r_syn(k) results.r_syn(k)], [0 yl(2)], 'k-', 'LineWidth',2);
    title(sprintf('%s vs %s (entry)', results.RegionA(k), results.RegionB(k)));
    xlabel('r'); ylabel('prob'); hold off;
end
saveas(figS, fullfile(out_dir, sprintf('synergy_nulls_%s_%s.png', activity_tag, condition)));
print(figS, fullfile(out_dir, sprintf('synergy_nulls_%s_%s.svg', activity_tag, condition)), '-dsvg');
close(figS);

fprintf('Done. Results saved to %s (CSV + figures).\n', outcsv);
end

function close_waitbar_safe(wb)
if isempty(wb)
    return;
end
try
    if isvalid(wb)
        close(wb);
    end
catch
end
end

%% Local helper functions
function [r, pperm, r_null] = perm_test_circshift(A, B, mask, Nperm)
if nargin < 4 || isempty(Nperm); Nperm = 1000; end
if nargin < 3 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
vecA = A; vecB = B;
r = matrix_corr(vecA, vecB, mask);
n = size(A,1);
r_null = nan(Nperm,1);
for t = 1:Nperm
    sh = randi(n);
    Ashift = circshift(vecA, sh, 1);
    Ashift = circshift(Ashift, sh, 2);
    r_null(t) = matrix_corr(Ashift, vecB, mask);
end
pperm = (sum(r_null >= r) + 1) / (Nperm + 1);
end

function [r, pperm, r_null] = perm_test_entry_shuffle(A,B,mask,Nperm)
if nargin < 5 || isempty(Nperm); Nperm = 1000; end
if nargin < 3 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
vecA = A(mask); vecB = B(mask);
valid = ~(isnan(vecA) | isnan(vecB)); vecA = vecA(valid); vecB = vecB(valid);
r = corr(vecA(:), vecB(:), 'Rows','complete');
n = numel(vecA);
r_null = nan(Nperm,1);
for t = 1:Nperm
    idx = randperm(n);
    r_null(t) = corr(vecA(idx), vecB, 'Rows','complete');
end
pperm = (sum(r_null >= r) + 1) / (Nperm + 1);
end

function [r_mean, CI] = bootstrap_similarity(A_all, B_all, mask, Nboot)
% A_all, B_all: time x time x nitems (nitems may be 0 -> return NaN)
if nargin < 4 || isempty(Nboot); Nboot = 1000; end
if isempty(A_all) || isempty(B_all)
    r_mean = NaN; CI = [NaN NaN]; return;
end
nA = size(A_all,3); nB = size(B_all,3);
% If no pairs, return NaN
if nA == 0 || nB == 0
    r_mean = NaN; CI = [NaN NaN]; return;
end
% bootstrap by resampling indices (resample smaller of the two if different sizes)
n = min(nA, nB);
r_boot = nan(Nboot,1);
for b = 1:Nboot
    ia = randi(nA, n, 1);
    ib = randi(nB, n, 1);
    Amean = nanmean(A_all(:,:,ia),3);
    Bmean = nanmean(B_all(:,:,ib),3);
    r_boot(b) = matrix_corr(Amean, Bmean, mask);
end
r_mean = mean(r_boot);
CI = prctile(r_boot, [2.5 97.5]);
end

function [r, p] = matrix_corr(A,B,mask)
if nargin < 3 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
vecA = A(mask); vecB = B(mask);
valid = ~(isnan(vecA) | isnan(vecB)); vecA = vecA(valid); vecB = vecB(valid);
if isempty(vecA)
    r = NaN; p = NaN; return;
end
[r,p] = corr(vecA(:), vecB(:), 'Type', 'Pearson', 'Rows', 'complete');
end

function s = matrix_ssim(A,B,mask)
if nargin < 3 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
A2 = A; B2 = B; A2(~mask) = 0; B2(~mask) = 0;
v = [A2(mask); B2(mask)];
minv = min(v); maxv = max(v);
if maxv > minv
    An = (A2 - minv) / (maxv - minv);
    Bn = (B2 - minv) / (maxv - minv);
else
    An = zeros(size(A2)); Bn = zeros(size(B2));
end
% Use MATLAB's SSIM when available; otherwise fall back to Pearson on
% the normalized matrices as a weaker proxy.
if exist('ssim','file') == 2
    s = ssim(An, Bn);
else
    s = matrix_corr(An, Bn, mask);
end
end

function [p, obs, r_null] = perm_test_circshift_metric(metricFcn, A, B, mask, Nperm)
% Generic circular-shift permutation test using a metric function handle
if nargin < 5 || isempty(Nperm); Nperm = 1000; end
if nargin < 4 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
obs = metricFcn(A, B, mask);
n = size(A,1);
r_null = nan(Nperm,1);
for t = 1:Nperm
    sh = randi(n);
    Ashift = circshift(A, sh, 1);
    Ashift = circshift(Ashift, sh, 2);
    r_null(t) = metricFcn(Ashift, B, mask);
end
p = (sum(r_null >= obs) + 1) / (Nperm + 1);
end

function [p, obs, r_null] = perm_test_entry_shuffle_metric(metricFcn, A, B, mask, Nperm)
% Generic entry-shuffle permutation test using a metric function handle
if nargin < 5 || isempty(Nperm); Nperm = 1000; end
if nargin < 4 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
vecA = A(mask); vecB = B(mask);
valid = ~(isnan(vecA) | isnan(vecB)); vecA = vecA(valid); vecB = vecB(valid);
obs = metricFcn(A, B, mask);
n = numel(vecA);
r_null = nan(Nperm,1);
for t = 1:Nperm
    idx = randperm(n);
    % build permuted matrix: copy A and replace masked entries with permuted values
    Aprime = A;
    tmp = A(mask);
    tmp(valid) = vecA(idx);
    Aprime(mask) = tmp;
    r_null(t) = metricFcn(Aprime, B, mask);
end
p = (sum(r_null >= obs) + 1) / (Nperm + 1);
end

function [ssim_mean, CI] = bootstrap_ssim(A_all, B_all, mask, Nboot)
% Bootstrap CI for SSIM computed from resampled mean matrices
if nargin < 4 || isempty(Nboot); Nboot = 1000; end
if isempty(A_all) || isempty(B_all)
    ssim_mean = NaN; CI = [NaN NaN]; return;
end
nA = size(A_all,3); nB = size(B_all,3);
if nA == 0 || nB == 0
    ssim_mean = NaN; CI = [NaN NaN]; return;
end
n = min(nA, nB);
ssim_boot = nan(Nboot,1);
for b = 1:Nboot
    ia = randi(nA, n, 1);
    ib = randi(nB, n, 1);
    Amean = nanmean(A_all(:,:,ia),3);
    Bmean = nanmean(B_all(:,:,ib),3);
    ssim_boot(b) = matrix_ssim(Amean, Bmean, mask);
end
ssim_mean = mean(ssim_boot);
CI = prctile(ssim_boot, [2.5 97.5]);
end
