function results = compare_drosophila_marmoset_coI(participants, basefold, datatype, activity_tag, condition, cutoff, timing, marmoset_csv_files, out_dir, scale_factor, marmoset_stim_onset_ms, marmoset_downsample_factor, Nperm, Nboot)
%COMPARE_DROSOPHILA_MARMOSET_COI Compare drosophila and scaled marmoset CoI matrices
%
% This function compares Drosophila CoI matrices (group-level averages across regions)
% with marmoset CoI matrices loaded from CSV files. Both are scaled to comparable
% ranges, then compared using Pearson correlation and SSIM metrics with permutation
% significance testing.
%
% INPUTS:
%   participants            - Cell array of fly names to include
%   basefold                - Base folder for data
%   datatype                - Type of data being analyzed
%   activity_tag            - Activity condition tag (e.g., 'wake_sleep')
%   condition               - Experimental condition (e.g., 'BSLEEP')
%   cutoff                  - Electrode cutoff for central/peripheral division
%   timing                  - Time vector for drosophila analysis
%   marmoset_csv_files      - Struct mapping region names to CSV file paths.
%                             Example: marmoset_csv_files.peripheral = 'path/to/peripheral.csv'
%                                      marmoset_csv_files.central = 'path/to/central.csv'
%                             Supported regions: peripheral, central, central_peripheral, 
%                             peripheral_peripheral, central_central
%   out_dir                 - Directory where outputs will be saved (default: Results/Cross_Species)
%   scale_factor            - Upscaling factor for Drosophila CoI matrix size (default: 3.0)
%   marmoset_stim_onset_ms  - Stim onset index for marmoset matrix in ms/index units (default: 100)
%   marmoset_downsample_factor - Factor by which marmoset CSV matrices were
%                             downsampled relative to original ms/index units
%                             (default: 1, no downsampling). This can be a
%                             scalar used for every region, or a struct with
%                             region-specific fields, e.g.
%                             struct('peripheral',1,'central',3).
%   Nperm                  - Number of permutations for circular-shift and
%                             entry-shuffle null tests (default: 1000)
%   Nboot                  - Number of bootstrap samples for confidence
%                             intervals (default: 1000)
%
% OUTPUT:
%   results                 - Structure containing comparison metrics, p-values, and CIs
%
% USAGE:
%   csv_map = struct('peripheral', 'path/to/Figure_3a_peripheral.csv', ...
%                    'central', 'path/to/Figure_3c_central.csv');
%   results = compare_drosophila_marmoset_coI(...
%       fly_list, basefold, 'LFP', 'wake_sleep', 'BSLEEP', 40, timing, ...
%       csv_map, output_dir, 3.0, 100);
%
% The function:
% 1. Aggregates Drosophila CoI data across flies and regions (like compare_region_patterns)
% 2. Loads marmoset CoI matrices from CSV files (one file per region)
% 3. Scales both to zero-mean unit variance for comparison
% 4. Computes Pearson r and SSIM metrics
% 5. Runs permutation tests (circular shift, entry shuffle) for significance
% 6. Bootstraps confidence intervals
% 7. Generates comparison figures and CSV summary

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fileparts(this_dir));

if nargin < 9 || isempty(out_dir)
    out_dir = fullfile(basefold, 'Results', 'Cross_Species_Comparison');
end
if nargin < 10 || isempty(scale_factor)
    scale_factor = 3.0;
end
if nargin < 11 || isempty(marmoset_stim_onset_ms)
    marmoset_stim_onset_ms = 100;
end
if nargin < 12 || isempty(marmoset_downsample_factor)
    marmoset_downsample_factor = 1;
end
if nargin < 13 || isempty(Nperm)
    Nperm = 1000;
end
if nargin < 14 || isempty(Nboot)
    Nboot = 1000;
end
if isnumeric(marmoset_downsample_factor) && any(marmoset_downsample_factor(:) <= 0)
    error('marmoset_downsample_factor must be positive.');
end
if ~isscalar(Nperm) || Nperm < 1 || Nperm ~= round(Nperm)
    error('Nperm must be a positive integer scalar.');
end
if ~isscalar(Nboot) || Nboot < 1 || Nboot ~= round(Nboot)
    error('Nboot must be a positive integer scalar.');
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Keep only post-0ms time points so comparison starts at (0,0) ms.
post0_idx = find(timing >= 0);
if isempty(post0_idx)
    error('No time points >= 0 ms found in timing vector.');
end
n_post0 = numel(post0_idx);
n_pre = sum(timing < 0);
n_all = numel(timing);
fprintf('Using post-0ms window: %d of %d time points (timing >= 0).\n', numel(post0_idx), numel(timing));
fprintf('Upscaling Drosophila CoI matrices by factor %.3f with linear interpolation.\n', scale_factor);
fprintf('Using marmoset stim onset index: %d.\n', round(marmoset_stim_onset_ms));
if isnumeric(marmoset_downsample_factor)
    fprintf('Using marmoset downsample factor: %.3f.\n', marmoset_downsample_factor);
else
    fprintf('Using region-specific marmoset downsample factors.\n');
end
fprintf('Using %d permutations and %d bootstrap samples.\n', Nperm, Nboot);

target_n = max(1, round(n_post0 * scale_factor));
target_n_all = max(1, round(n_all * scale_factor));
target_pre = max(0, round(n_pre * scale_factor));
fprintf('Comparison window size: %d x %d (post-stim only).\n', target_n, target_n);
fprintf('Comparison window size with pre-stim: %d x %d (includes %d pre-stim samples after scaling).\n', target_n_all, target_n_all, target_pre);

%regions = {'peripheral', 'central', 'central_peripheral', 'peripheral_peripheral', 'central_central'};
%region_labels = {'Peripheral', 'Central', 'Central-Peripheral', 'Peripheral-Peripheral', 'Central-Central'};
regions = {'peripheral', 'central'};
region_labels = {'Peripheral', 'Central'};
n_regions = numel(regions);

% -----------------------------------------------------------------
% STEP 1: Aggregate Drosophila CoI data across flies (same as compare_region_patterns)
% -----------------------------------------------------------------
fprintf('Aggregating Drosophila CoI data across %d flies...\n', numel(participants));

group = struct();
for r = 1:n_regions
    group.(regions{r}).FFi_all = [];
    group.(regions{r}).mask_all = [];
    group.(regions{r}).mi_all = [];
    group.(regions{r}).n_pairs = 0;
end

% Loop over participants and aggregate
for s = 1:numel(participants)
    participant_name = char(participants(s));
    flyData = max_get_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, timing);
    for r = 1:n_regions
        reg = regions{r};
        fly = flyData.(reg);
        if isempty(fly.FFi_all)
            continue;
        end
        group.(reg).FFi_all = cat(3, group.(reg).FFi_all, fly.FFi_all);
        group.(reg).mask_all = cat(3, group.(reg).mask_all, fly.mask);
        if isfield(fly, 'mi1') && ~isempty(fly.mi1)
            group.(reg).mi_all = [group.(reg).mi_all, fly.mi1];
        end
        if isfield(fly, 'mi2') && ~isempty(fly.mi2)
            group.(reg).mi_all = [group.(reg).mi_all, fly.mi2];
        end
        group.(reg).n_pairs = group.(reg).n_pairs + size(fly.FFi_all, 3);
    end
end

% Compute group means and masks
dros_struct = struct();
dros_struct_prestim = struct();
for r = 1:n_regions
    reg = regions{r};
    if isempty(group.(reg).FFi_all)
        warning('Region %s has no data', reg);
        dros_struct(r).mean = nan(target_n, target_n);
        dros_struct(r).mask = false(target_n, target_n);
        dros_struct(r).mi = nan(target_n, 1);
        dros_struct(r).mi_sem = nan(target_n, 1);
        dros_struct_prestim(r).mean = nan(target_n_all, target_n_all);
        dros_struct_prestim(r).mask = false(target_n_all, target_n_all);
        dros_struct_prestim(r).mi = nan(target_n_all, 1);
        dros_struct_prestim(r).mi_sem = nan(target_n_all, 1);
        dros_struct(r).n_pairs = 0;
        dros_struct_prestim(r).n_pairs = 0;
    else
        dros_mean_all = nanmean(group.(reg).FFi_all, 3);
        dros_mask_all = mean(group.(reg).mask_all, 3, 'omitnan') > 0;

        dros_mean = dros_mean_all;
        dros_mask = dros_mask_all;

        % Restrict to post-0ms in both dimensions.
        dros_mean = dros_mean(post0_idx, post0_idx);
        dros_mask = dros_mask(post0_idx, post0_idx);

        % Upscale Drosophila matrix and linearly interpolate missing values.
        [dros_mean_up, dros_mask_up] = upscale_and_interpolate_matrix(dros_mean, dros_mask, scale_factor);

        dros_struct(r).mean = dros_mean_up;
        dros_struct(r).mask = dros_mask_up;
        [dros_struct(r).mi, dros_struct(r).mi_sem] = build_region_mi_curve(group.(reg).mi_all, post0_idx, target_n);
        dros_struct(r).n_pairs = group.(reg).n_pairs;

        % Also build pre+post matrix for additional comparison.
        [dros_mean_up_all, dros_mask_up_all] = upscale_and_interpolate_matrix(dros_mean_all, dros_mask_all, scale_factor);
        dros_struct_prestim(r).mean = dros_mean_up_all;
        dros_struct_prestim(r).mask = dros_mask_up_all;
        [dros_struct_prestim(r).mi, dros_struct_prestim(r).mi_sem] = build_region_mi_curve(group.(reg).mi_all, 1:n_all, target_n_all);
        dros_struct_prestim(r).n_pairs = group.(reg).n_pairs;
    end
    dros_struct(r).name = reg;
    dros_struct(r).label = region_labels{r};
    dros_struct_prestim(r).name = reg;
    dros_struct_prestim(r).label = region_labels{r};
end

% -----------------------------------------------------------------
% STEP 2: Load marmoset CoI data from CSV files (one per region)
% -----------------------------------------------------------------
fprintf('Loading marmoset CoI data from CSV files...\n');

marm_struct = struct();
marm_struct_prestim = struct();
for r = 1:n_regions
    reg = regions{r};
    
    % Check if this region has a CSV file in the mapping
    if isstruct(marmoset_csv_files) && isfield(marmoset_csv_files, reg)
        csv_file = marmoset_csv_files.(reg);
    else
        csv_file = [];
    end
    
    if isempty(csv_file)
        fprintf('  %s: No CSV file specified, skipping...\n', reg);
        marm_struct(r).mean = nan(size(dros_struct(r).mean));
        marm_struct(r).mask = false(size(dros_struct(r).mean));
        marm_struct(r).label = region_labels{r};
        marm_struct(r).name = reg;
        marm_struct_prestim(r).mean = nan(size(dros_struct_prestim(r).mean));
        marm_struct_prestim(r).mask = false(size(dros_struct_prestim(r).mean));
        marm_struct_prestim(r).label = region_labels{r};
        marm_struct_prestim(r).name = reg;
        continue;
    end
    
    try
        % Load from CSV file
        marm_array = readmatrix(csv_file);
        region_downsample_factor = get_region_downsample_factor(marmoset_downsample_factor, reg);

        % Keep only marmoset post-stim matrix from the provided onset.
        % If the CSV is downsampled, convert original ms/index units into
        % saved matrix indices, crop in the saved matrix, then interpolate
        % back to the Drosophila comparison grid.
        onset_idx = max(1, round(marmoset_stim_onset_ms / region_downsample_factor));
        target_n_marm = max(1, round(target_n / region_downsample_factor));
        target_n_all_marm = max(1, round(target_n_all / region_downsample_factor));
        target_pre_marm = max(0, round(target_pre / region_downsample_factor));
        if size(marm_array, 1) < onset_idx || size(marm_array, 2) < onset_idx
            error('Marmoset matrix is smaller than onset index %d (%d x %d).', onset_idx, size(marm_array,1), size(marm_array,2));
        end
        marm_post = marm_array(onset_idx:end, onset_idx:end);

        % Compare only the fixed post-stim window matching scaled Drosophila extent.
        if size(marm_post, 1) < target_n_marm || size(marm_post, 2) < target_n_marm
            error('Post-onset marmoset matrix too small for %d x %d saved-sample window (available: %d x %d).', ...
                target_n_marm, target_n_marm, size(marm_post,1), size(marm_post,2));
        end
        marm_crop_native = marm_post(1:target_n_marm, 1:target_n_marm);
        marm_crop = interp2_resize_linear(marm_crop_native, [target_n, target_n]);

        % Additional window including pre-stim time.
        pre_start_idx = max(1, onset_idx - target_pre_marm);
        effective_pre = onset_idx - pre_start_idx;
        if effective_pre < target_pre_marm
            fprintf('  Note [%s]: requested pre-stim=%d saved samples but available=%d at onset=%d (start clamped to 1).\n', ...
                reg, target_pre_marm, effective_pre, onset_idx);
        end
        end_idx_all = pre_start_idx + target_n_all_marm - 1;
        if end_idx_all > size(marm_array, 1) || end_idx_all > size(marm_array, 2)
            error('Marmoset matrix too small for pre+post saved-sample window %d x %d from start index %d (matrix: %d x %d).', ...
                target_n_all_marm, target_n_all_marm, pre_start_idx, size(marm_array,1), size(marm_array,2));
        end
        marm_crop_prestim_native = marm_array(pre_start_idx:end_idx_all, pre_start_idx:end_idx_all);
        marm_crop_prestim = interp2_resize_linear(marm_crop_prestim_native, [target_n_all, target_n_all]);

        marm_struct(r).mean = marm_crop;
        marm_struct(r).mask = ~isnan(marm_crop);
        marm_struct(r).label = region_labels{r};
        marm_struct(r).name = reg;
        marm_struct_prestim(r).mean = marm_crop_prestim;
        marm_struct_prestim(r).mask = ~isnan(marm_crop_prestim);
        marm_struct_prestim(r).label = region_labels{r};
        marm_struct_prestim(r).name = reg;
        fprintf('  %s <- %s: using post-onset %d x %d saved-sample window, resized to %d x %d\n', ...
            reg, csv_file, size(marm_crop_native, 1), size(marm_crop_native, 2), size(marm_crop, 1), size(marm_crop, 2));
        fprintf('    Sanity [%s]: dros_post0=%dx%d, dros_scaled=%dx%d, marm_full=%dx%d, marm_post_onset=%dx%d, compared=%dx%d, downsample=%.3f\n', ...
            reg, n_post0, n_post0, size(dros_struct(r).mean,1), size(dros_struct(r).mean,2), ...
            size(marm_array,1), size(marm_array,2), size(marm_post,1), size(marm_post,2), ...
            size(marm_crop,1), size(marm_crop,2), region_downsample_factor);
        fprintf('    Sanity+Pre [%s]: dros_all_scaled=%dx%d, marm_start_idx=%d, native_with_pre=%dx%d, compared_with_pre=%dx%d\n', ...
            reg, size(dros_struct_prestim(r).mean,1), size(dros_struct_prestim(r).mean,2), pre_start_idx, ...
            size(marm_crop_prestim_native,1), size(marm_crop_prestim_native,2), ...
            size(marm_crop_prestim,1), size(marm_crop_prestim,2));
    catch ME
        fprintf('  Warning: Could not load %s from CSV %s: %s\n', reg, csv_file, ME.message);
        marm_struct(r).mean = nan(size(dros_struct(r).mean));
        marm_struct(r).mask = false(size(dros_struct(r).mean));
        marm_struct(r).label = region_labels{r};
        marm_struct(r).name = reg;
        marm_struct_prestim(r).mean = nan(size(dros_struct_prestim(r).mean));
        marm_struct_prestim(r).mask = false(size(dros_struct_prestim(r).mean));
        marm_struct_prestim(r).label = region_labels{r};
        marm_struct_prestim(r).name = reg;
        fprintf('    Sanity [%s]: dros_post0=%dx%d, dros_scaled=%dx%d, marm_full=NA, marm_post_onset=NA, compared=%dx%d\n', ...
            reg, n_post0, n_post0, size(dros_struct(r).mean,1), size(dros_struct(r).mean,2), ...
            size(dros_struct(r).mean,1), size(dros_struct(r).mean,2));
        fprintf('    Sanity+Pre [%s]: dros_all_scaled=%dx%d, marm_start_idx=NA, compared_with_pre=%dx%d\n', ...
            reg, size(dros_struct_prestim(r).mean,1), size(dros_struct_prestim(r).mean,2), ...
            size(dros_struct_prestim(r).mean,1), size(dros_struct_prestim(r).mean,2));
    end
end

% -----------------------------------------------------------------
% STEP 3: Scale both datasets (zero-mean, unit variance)
% -----------------------------------------------------------------
fprintf('Scaling Drosophila and marmoset matrices...\n');

for r = 1:n_regions
    % Drosophila
    if ~all(isnan(dros_struct(r).mean(:)))
        dros_vec = dros_struct(r).mean(dros_struct(r).mask);
        dros_vec = dros_vec(~isnan(dros_vec));
        if ~isempty(dros_vec)
            dros_scaled = (dros_struct(r).mean - nanmean(dros_vec)) / nanstd(dros_vec);
        else
            dros_scaled = dros_struct(r).mean;
        end
    else
        dros_scaled = dros_struct(r).mean;
    end
    dros_struct(r).scaled = dros_scaled;

    % Drosophila (including pre-stim)
    if ~all(isnan(dros_struct_prestim(r).mean(:)))
        dros_vec_all = dros_struct_prestim(r).mean(dros_struct_prestim(r).mask);
        dros_vec_all = dros_vec_all(~isnan(dros_vec_all));
        if ~isempty(dros_vec_all)
            dros_scaled_all = (dros_struct_prestim(r).mean - nanmean(dros_vec_all)) / nanstd(dros_vec_all);
        else
            dros_scaled_all = dros_struct_prestim(r).mean;
        end
    else
        dros_scaled_all = dros_struct_prestim(r).mean;
    end
    dros_struct_prestim(r).scaled = dros_scaled_all;
    
    % Marmoset
    if ~all(isnan(marm_struct(r).mean(:)))
        marm_vec = marm_struct(r).mean(marm_struct(r).mask);
        marm_vec = marm_vec(~isnan(marm_vec));
        if ~isempty(marm_vec)
            marm_scaled = (marm_struct(r).mean - nanmean(marm_vec)) / nanstd(marm_vec);
        else
            marm_scaled = marm_struct(r).mean;
        end
    else
        marm_scaled = marm_struct(r).mean;
    end
    marm_struct(r).scaled = marm_scaled;

    % Marmoset (including pre-stim)
    if ~all(isnan(marm_struct_prestim(r).mean(:)))
        marm_vec_all = marm_struct_prestim(r).mean(marm_struct_prestim(r).mask);
        marm_vec_all = marm_vec_all(~isnan(marm_vec_all));
        if ~isempty(marm_vec_all)
            marm_scaled_all = (marm_struct_prestim(r).mean - nanmean(marm_vec_all)) / nanstd(marm_vec_all);
        else
            marm_scaled_all = marm_struct_prestim(r).mean;
        end
    else
        marm_scaled_all = marm_struct_prestim(r).mean;
    end
    marm_struct_prestim(r).scaled = marm_scaled_all;
end

% -----------------------------------------------------------------
% STEP 3.1: Load/align marginal MI timecourses for matrix panels
% -----------------------------------------------------------------
dros_time_post = linspace(timing(post0_idx(1)), timing(post0_idx(end)), target_n);
dros_time_prestim = linspace(timing(1), timing(end), target_n_all);
marmMI = load_marmoset_mi_curves(basefold, target_n, target_n_all, target_pre, marmoset_stim_onset_ms);
marm_time_post = marmMI.time_post;
marm_time_prestim = marmMI.time_with_prestim;
for r = 1:n_regions
    marm_struct(r).mi = marmMI.(regions{r}).post;
    marm_struct(r).mi_sem = marmMI.(regions{r}).post_sem;
    marm_struct_prestim(r).mi = marmMI.(regions{r}).with_prestim;
    marm_struct_prestim(r).mi_sem = marmMI.(regions{r}).with_prestim_sem;
end

% -----------------------------------------------------------------
% STEP 3.5: QC plot for normalized cropped Drosophila matrices
% -----------------------------------------------------------------
fprintf('Generating Drosophila QC plot (normalized, cropped/scaled)...\n');

all_dros_scaled = [];
for r = 1:n_regions
    if ~isempty(dros_struct(r).scaled)
        all_dros_scaled = [all_dros_scaled; dros_struct(r).scaled(:)];
    end
end
all_dros_scaled = all_dros_scaled(isfinite(all_dros_scaled));

if isempty(all_dros_scaled)
    dros_scaled_clim = [-2 2];
else
    dros_scaled_absmax = max(abs([min(all_dros_scaled), max(all_dros_scaled)]));
    if dros_scaled_absmax == 0
        dros_scaled_absmax = 1;
    end
    dros_scaled_clim = [-dros_scaled_absmax dros_scaled_absmax];
end

fig_dros_qc = figure('Position', [100 100 1400 500], 'Renderer', 'painters');
sgtitle(sprintf('Drosophila CoI QC: Normalized, Cropped, and Scaled (%s, %s)', activity_tag, condition));
tiledlayout(1, n_regions, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for r = 1:n_regions
    nexttile(r);
    imagesc(dros_struct(r).scaled);
    set(gca, 'YDir', 'normal');
    axis square;
    colorbar;
    colormap(gca, redblue(256));
    caxis(dros_scaled_clim);
    title(sprintf('Normalized - %s', region_labels{r}));
end

saveas(fig_dros_qc, fullfile(out_dir, sprintf('drosophila_qc_normalized_cropped_scaled_%s_%s.png', activity_tag, condition)));
print_editable_svg(fig_dros_qc, fullfile(out_dir, sprintf('drosophila_qc_normalized_cropped_scaled_%s_%s.svg', activity_tag, condition)));
close(fig_dros_qc);

% -----------------------------------------------------------------
% STEP 3.6: QC plot for marmoset loading (raw vs normalized)
% -----------------------------------------------------------------
fprintf('Generating marmoset QC plot (raw and normalized)...\n');

all_marm_raw = [];
all_marm_scaled = [];
for r = 1:n_regions
    if ~isempty(marm_struct(r).mean)
        all_marm_raw = [all_marm_raw; marm_struct(r).mean(:)];
        all_marm_scaled = [all_marm_scaled; marm_struct(r).scaled(:)];
    end
end
all_marm_raw = all_marm_raw(isfinite(all_marm_raw));
all_marm_scaled = all_marm_scaled(isfinite(all_marm_scaled));

if isempty(all_marm_raw)
    raw_clim = [-1 1];
else
    raw_absmax = max(abs([min(all_marm_raw), max(all_marm_raw)]));
    if raw_absmax == 0
        raw_absmax = 1;
    end
    raw_clim = [-raw_absmax raw_absmax];
end

if isempty(all_marm_scaled)
    scaled_clim = [-2 2];
else
    scaled_absmax = max(abs([min(all_marm_scaled), max(all_marm_scaled)]));
    if scaled_absmax == 0
        scaled_absmax = 1;
    end
    scaled_clim = [-scaled_absmax scaled_absmax];
end

fig_qc = figure('Position', [100 100 1600 700], 'Renderer', 'painters');
sgtitle(sprintf('Marmoset CoI QC: Raw vs Normalized (%s, %s)', activity_tag, condition));
tiledlayout(2, n_regions, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for r = 1:n_regions
    % Top row: raw marmoset matrices
    nexttile(r);
    imagesc(marm_struct(r).mean);
    set(gca, 'YDir', 'normal');
    axis square;
    colorbar;
    colormap(gca, redblue(256));
    caxis(raw_clim);
    title(sprintf('Raw - %s', region_labels{r}));

    % Bottom row: normalized marmoset matrices
    nexttile(n_regions + r);
    imagesc(marm_struct(r).scaled);
    set(gca, 'YDir', 'normal');
    axis square;
    colorbar;
    colormap(gca, redblue(256));
    caxis(scaled_clim);
    title(sprintf('Normalized - %s', region_labels{r}));
end

saveas(fig_qc, fullfile(out_dir, sprintf('marmoset_qc_raw_vs_normalized_%s_%s.png', activity_tag, condition)));
print_editable_svg(fig_qc, fullfile(out_dir, sprintf('marmoset_qc_raw_vs_normalized_%s_%s.svg', activity_tag, condition)));
close(fig_qc);

% -----------------------------------------------------------------
% STEP 4: Compare metrics per region (Pearson r, SSIM)
% -----------------------------------------------------------------
fprintf('Computing comparison metrics...\n');

% Create output table
results = table;
results.Region = strings(n_regions, 1);
results.Label = strings(n_regions, 1);
results.N_dros_pairs = nan(n_regions, 1);
results.r_pearson = nan(n_regions, 1);
results.p_circ = nan(n_regions, 1);
results.p_entry = nan(n_regions, 1);
results.CI_low = nan(n_regions, 1);
results.CI_high = nan(n_regions, 1);
results.r_redundancy = nan(n_regions, 1);
results.p_circ_redundancy = nan(n_regions, 1);
results.p_entry_redundancy = nan(n_regions, 1);
results.CI_low_redundancy = nan(n_regions, 1);
results.CI_high_redundancy = nan(n_regions, 1);
results.r_synergy = nan(n_regions, 1);
results.p_circ_synergy = nan(n_regions, 1);
results.p_entry_synergy = nan(n_regions, 1);
results.CI_low_synergy = nan(n_regions, 1);
results.CI_high_synergy = nan(n_regions, 1);
results.ssim = nan(n_regions, 1);
results.ssim_p_circ = nan(n_regions, 1);
results.ssim_q_circ = nan(n_regions, 1);
results.ssim_p_entry = nan(n_regions, 1);
results.ssim_q_entry = nan(n_regions, 1);
results.ssim_CI_low = nan(n_regions, 1);
results.ssim_CI_high = nan(n_regions, 1);

% Storage for null distributions (for plotting)
nulls = struct();
ssim_debug = struct();
for r = 1:n_regions
    nulls(r).circ_r = nan(Nperm, 1);
    nulls(r).entry_r = nan(Nperm, 1);
    nulls(r).circ_ssim = nan(Nperm, 1);
    nulls(r).entry_ssim = nan(Nperm, 1);
    ssim_debug(r).region = regions{r};
    ssim_debug(r).X = [];
    ssim_debug(r).Y = [];
    ssim_debug(r).mask = [];
    ssim_debug(r).ssim = NaN;
end

for r = 1:n_regions
    results.Region(r) = dros_struct(r).name;
    results.Label(r) = dros_struct(r).label;
    results.N_dros_pairs(r) = dros_struct(r).n_pairs;
    
    % Get scaled matrices
    X = dros_struct(r).scaled;
    Y = marm_struct(r).scaled;

    if ~isequal(size(X), size(Y))
        error('Size mismatch in region %s after windowing: drosophila=%dx%d, marmoset=%dx%d', ...
            results.Region(r), size(X,1), size(X,2), size(Y,1), size(Y,2));
    end
    
    % Union mask (non-NaN in both)
    mask = ~isnan(X) & ~isnan(Y);
    if all(~mask(:))
        % Fallback to all finite cells
        mask = isfinite(X) & isfinite(Y);
    end
    
    if ~any(mask(:))
        fprintf('  Skipping region %s: no valid data\n', results.Region(r));
        continue;
    end
    
    % Pearson correlation
    [r_obs, p_circ, r_null_circ] = perm_test_circshift(X, Y, mask, Nperm);
    [~, p_entry, r_null_entry] = perm_test_entry_shuffle(X, Y, mask, Nperm);
    [~, CI] = bootstrap_similarity_scalar(X, Y, mask, Nboot);
    
    results.r_pearson(r) = r_obs;
    results.p_circ(r) = p_circ;
    results.p_entry(r) = p_entry;
    results.CI_low(r) = CI(1);
    results.CI_high(r) = CI(2);

    % Pearson correlation split into redundancy and synergy components.
    % These use the matched raw CoI matrices, not the z-scored signed matrices:
    % redundancy = positive CoI; synergy = absolute magnitude of negative CoI.
    component_metrics = compute_component_pearson_metrics(dros_struct(r).mean, marm_struct(r).mean, Nperm, Nboot);
    results.r_redundancy(r) = component_metrics.redundancy.r;
    results.p_circ_redundancy(r) = component_metrics.redundancy.p_circ;
    results.p_entry_redundancy(r) = component_metrics.redundancy.p_entry;
    results.CI_low_redundancy(r) = component_metrics.redundancy.CI(1);
    results.CI_high_redundancy(r) = component_metrics.redundancy.CI(2);
    results.r_synergy(r) = component_metrics.synergy.r;
    results.p_circ_synergy(r) = component_metrics.synergy.p_circ;
    results.p_entry_synergy(r) = component_metrics.synergy.p_entry;
    results.CI_low_synergy(r) = component_metrics.synergy.CI(1);
    results.CI_high_synergy(r) = component_metrics.synergy.CI(2);
    
    % SSIM
    ssim_obs = matrix_ssim(X, Y, mask);
    [ssim_p_circ, ~, ssim_null_circ] = perm_test_circshift_metric(@matrix_ssim, X, Y, mask, Nperm);
    [ssim_p_entry, ~, ssim_null_entry] = perm_test_entry_shuffle_metric(@matrix_ssim, X, Y, mask, Nperm);
    [~, ssim_CI] = bootstrap_ssim_scalar(X, Y, mask, Nboot);
    
    results.ssim(r) = ssim_obs;
    results.ssim_p_circ(r) = ssim_p_circ;
    results.ssim_p_entry(r) = ssim_p_entry;
    results.ssim_CI_low(r) = ssim_CI(1);
    results.ssim_CI_high(r) = ssim_CI(2);
    ssim_debug(r).X = X;
    ssim_debug(r).Y = Y;
    ssim_debug(r).mask = mask;
    ssim_debug(r).ssim = ssim_obs;
    
    % Store null distributions
    nulls(r).circ_r = r_null_circ;
    nulls(r).entry_r = r_null_entry;
    nulls(r).circ_ssim = ssim_null_circ;
    nulls(r).entry_ssim = ssim_null_entry;
    
    fprintf('  %s: r=%.3f (p_circ=%.4f), r_red=%.3f, r_syn=%.3f, SSIM=%.3f (p_circ=%.4f)\n', ...
        results.Label(r), r_obs, p_circ, results.r_redundancy(r), results.r_synergy(r), ssim_obs, ssim_p_circ);
end

% -----------------------------------------------------------------
% STEP 4.5: Compare metrics per region INCLUDING pre-stim window
% -----------------------------------------------------------------
fprintf('Computing comparison metrics (including pre-stim window)...\n');

results_with_prestim = table;
results_with_prestim.Region = strings(n_regions, 1);
results_with_prestim.Label = strings(n_regions, 1);
results_with_prestim.N_dros_pairs = nan(n_regions, 1);
results_with_prestim.r_pearson = nan(n_regions, 1);
results_with_prestim.p_circ = nan(n_regions, 1);
results_with_prestim.p_entry = nan(n_regions, 1);
results_with_prestim.CI_low = nan(n_regions, 1);
results_with_prestim.CI_high = nan(n_regions, 1);
results_with_prestim.r_redundancy = nan(n_regions, 1);
results_with_prestim.p_circ_redundancy = nan(n_regions, 1);
results_with_prestim.p_entry_redundancy = nan(n_regions, 1);
results_with_prestim.CI_low_redundancy = nan(n_regions, 1);
results_with_prestim.CI_high_redundancy = nan(n_regions, 1);
results_with_prestim.r_synergy = nan(n_regions, 1);
results_with_prestim.p_circ_synergy = nan(n_regions, 1);
results_with_prestim.p_entry_synergy = nan(n_regions, 1);
results_with_prestim.CI_low_synergy = nan(n_regions, 1);
results_with_prestim.CI_high_synergy = nan(n_regions, 1);
results_with_prestim.ssim = nan(n_regions, 1);
results_with_prestim.ssim_p_circ = nan(n_regions, 1);
results_with_prestim.ssim_q_circ = nan(n_regions, 1);
results_with_prestim.ssim_p_entry = nan(n_regions, 1);
results_with_prestim.ssim_q_entry = nan(n_regions, 1);
results_with_prestim.ssim_CI_low = nan(n_regions, 1);
results_with_prestim.ssim_CI_high = nan(n_regions, 1);

for r = 1:n_regions
    results_with_prestim.Region(r) = dros_struct_prestim(r).name;
    results_with_prestim.Label(r) = dros_struct_prestim(r).label;
    results_with_prestim.N_dros_pairs(r) = dros_struct_prestim(r).n_pairs;

    X = dros_struct_prestim(r).scaled;
    Y = marm_struct_prestim(r).scaled;

    if ~isequal(size(X), size(Y))
        error('Size mismatch in region %s for pre-stim comparison: drosophila=%dx%d, marmoset=%dx%d', ...
            results_with_prestim.Region(r), size(X,1), size(X,2), size(Y,1), size(Y,2));
    end

    mask = ~isnan(X) & ~isnan(Y);
    if all(~mask(:))
        mask = isfinite(X) & isfinite(Y);
    end
    if ~any(mask(:))
        fprintf('  Skipping pre-stim region %s: no valid data\n', results_with_prestim.Region(r));
        continue;
    end

    [r_obs, p_circ, ~] = perm_test_circshift(X, Y, mask, Nperm);
    [~, p_entry, ~] = perm_test_entry_shuffle(X, Y, mask, Nperm);
    [~, CI] = bootstrap_similarity_scalar(X, Y, mask, Nboot);

    ssim_obs = matrix_ssim(X, Y, mask);
    [ssim_p_circ, ~, ~] = perm_test_circshift_metric(@matrix_ssim, X, Y, mask, Nperm);
    [ssim_p_entry, ~, ~] = perm_test_entry_shuffle_metric(@matrix_ssim, X, Y, mask, Nperm);
    [~, ssim_CI] = bootstrap_ssim_scalar(X, Y, mask, Nboot);

    results_with_prestim.r_pearson(r) = r_obs;
    results_with_prestim.p_circ(r) = p_circ;
    results_with_prestim.p_entry(r) = p_entry;
    results_with_prestim.CI_low(r) = CI(1);
    results_with_prestim.CI_high(r) = CI(2);

    component_metrics = compute_component_pearson_metrics(dros_struct_prestim(r).mean, marm_struct_prestim(r).mean, Nperm, Nboot);
    results_with_prestim.r_redundancy(r) = component_metrics.redundancy.r;
    results_with_prestim.p_circ_redundancy(r) = component_metrics.redundancy.p_circ;
    results_with_prestim.p_entry_redundancy(r) = component_metrics.redundancy.p_entry;
    results_with_prestim.CI_low_redundancy(r) = component_metrics.redundancy.CI(1);
    results_with_prestim.CI_high_redundancy(r) = component_metrics.redundancy.CI(2);
    results_with_prestim.r_synergy(r) = component_metrics.synergy.r;
    results_with_prestim.p_circ_synergy(r) = component_metrics.synergy.p_circ;
    results_with_prestim.p_entry_synergy(r) = component_metrics.synergy.p_entry;
    results_with_prestim.CI_low_synergy(r) = component_metrics.synergy.CI(1);
    results_with_prestim.CI_high_synergy(r) = component_metrics.synergy.CI(2);
    results_with_prestim.ssim(r) = ssim_obs;
    results_with_prestim.ssim_p_circ(r) = ssim_p_circ;
    results_with_prestim.ssim_p_entry(r) = ssim_p_entry;
    results_with_prestim.ssim_CI_low(r) = ssim_CI(1);
    results_with_prestim.ssim_CI_high(r) = ssim_CI(2);

    fprintf('  %s [with pre]: r=%.3f (p_circ=%.4f), r_red=%.3f, r_syn=%.3f, SSIM=%.3f (p_circ=%.4f)\n', ...
        results_with_prestim.Label(r), r_obs, p_circ, results_with_prestim.r_redundancy(r), results_with_prestim.r_synergy(r), ssim_obs, ssim_p_circ);
end

% FDR-correct SSIM p-values only across the pre-stim-included region tests.
% Post-only q-values are left as NaN because those analyses are no longer
% part of the planned inference family.
results_with_prestim.ssim_q_circ = bh_fdr(results_with_prestim.ssim_p_circ);
results_with_prestim.ssim_q_entry = bh_fdr(results_with_prestim.ssim_p_entry);

% -----------------------------------------------------------------
% STEP 4.6: SSIM interaction between region and cross-species similarity
% -----------------------------------------------------------------
fprintf('Computing SSIM region interaction: Peripheral minus Central...\n');
peripheral_result_idx = find(strcmp(results.Region, 'peripheral'), 1);
central_result_idx = find(strcmp(results.Region, 'central'), 1);
peripheral_result_idx_prestim = find(strcmp(results_with_prestim.Region, 'peripheral'), 1);
central_result_idx_prestim = find(strcmp(results_with_prestim.Region, 'central'), 1);

interaction_post = compute_ssim_region_interaction( ...
    dros_struct, marm_struct, 'post_only', Nperm, ...
    results.ssim(peripheral_result_idx), results.ssim(central_result_idx), ssim_debug);
interaction_with_prestim = compute_ssim_region_interaction( ...
    dros_struct_prestim, marm_struct_prestim, 'with_prestim', Nperm, ...
    results_with_prestim.ssim(peripheral_result_idx_prestim), results_with_prestim.ssim(central_result_idx_prestim));
interaction_results = [interaction_post.table; interaction_with_prestim.table];
interaction_results.QCircGreater = nan(height(interaction_results), 1);
interaction_results.QEntryGreater = nan(height(interaction_results), 1);
interaction_results.QCircTwoSided = nan(height(interaction_results), 1);
interaction_results.QEntryTwoSided = nan(height(interaction_results), 1);
with_prestim_interaction = interaction_results.ComparisonMode == "with_prestim";
interaction_results.QCircGreater(with_prestim_interaction) = bh_fdr(interaction_results.PCircGreater(with_prestim_interaction));
interaction_results.QEntryGreater(with_prestim_interaction) = bh_fdr(interaction_results.PEntryGreater(with_prestim_interaction));
interaction_results.QCircTwoSided(with_prestim_interaction) = bh_fdr(interaction_results.PCircTwoSided(with_prestim_interaction));
interaction_results.QEntryTwoSided(with_prestim_interaction) = bh_fdr(interaction_results.PEntryTwoSided(with_prestim_interaction));
plot_ssim_interaction_nulls(interaction_post, interaction_with_prestim, out_dir, activity_tag, condition);

function [A_up, mask_up] = upscale_and_interpolate_matrix(A, mask, scale_factor)
% Fill missing entries linearly and upscale to larger matrix dimensions.
    if nargin < 3 || isempty(scale_factor)
        scale_factor = 3.0;
    end
    if nargin < 2 || isempty(mask)
        mask = ~isnan(A);
    end

    A = double(A);
    mask = logical(mask);

    % Fill missing values using linear scattered interpolation on index grid.
    [rr, cc] = ndgrid(1:size(A,1), 1:size(A,2));
    valid = mask & isfinite(A);
    if any(valid(:))
        F_lin = scatteredInterpolant(rr(valid), cc(valid), A(valid), 'linear', 'nearest');
        A_filled = F_lin(rr, cc);
    else
        A_filled = zeros(size(A));
    end

    new_n1 = max(1, round(size(A_filled,1) * scale_factor));
    new_n2 = max(1, round(size(A_filled,2) * scale_factor));
    A_up = interp2_resize_linear(A_filled, [new_n1, new_n2]);

    mask_up = interp2_resize_linear(double(mask), [new_n1, new_n2]) > 0.1;
end

function B = interp2_resize_linear(A, new_size)
% Resize 2D matrix using linear interpolation on index coordinates.
    if isempty(A)
        B = nan(new_size);
        return;
    end

    n1 = size(A,1);
    n2 = size(A,2);
    r_new = linspace(1, n1, new_size(1));
    c_new = linspace(1, n2, new_size(2));
    [cc, rr] = meshgrid(1:n2, 1:n1);
    [cc_new, rr_new] = meshgrid(c_new, r_new);
    B = interp2(cc, rr, double(A), cc_new, rr_new, 'linear');
end

function downsample_factor = get_region_downsample_factor(downsample_setting, region_name)
% Return scalar downsample factor for a given marmoset region.
    if isnumeric(downsample_setting)
        downsample_factor = downsample_setting;
    elseif isstruct(downsample_setting) && isfield(downsample_setting, region_name)
        downsample_factor = downsample_setting.(region_name);
    elseif isstruct(downsample_setting)
        downsample_factor = 1;
        fprintf('  Note [%s]: no region-specific downsample factor found; using 1.\n', region_name);
    else
        error('marmoset_downsample_factor must be a positive scalar or a struct of positive region factors.');
    end

    if ~isscalar(downsample_factor) || ~isnumeric(downsample_factor) || downsample_factor <= 0
        error('Invalid marmoset downsample factor for region %s.', region_name);
    end
end

function component_metrics = compute_component_pearson_metrics(A, B, Nperm, Nboot)
% Pearson similarity for redundancy and synergy components separately.
%
% Redundancy is the positive CoI component. Synergy is converted to positive
% magnitude, i.e. -negative CoI, so the comparison asks whether the shape of
% the synergetic pattern is similar across species.
    A_red = A .* (A > 0);
    B_red = B .* (B > 0);
    A_syn = -A .* (A < 0);
    B_syn = -B .* (B < 0);

    component_metrics.redundancy = compute_one_component_pearson(A_red, B_red, Nperm, Nboot);
    component_metrics.synergy = compute_one_component_pearson(A_syn, B_syn, Nperm, Nboot);
end

function metric = compute_one_component_pearson(A, B, Nperm, Nboot)
    mask = isfinite(A) & isfinite(B) & ((A ~= 0) | (B ~= 0));
    if ~any(mask(:))
        mask = isfinite(A) & isfinite(B);
    end

    if ~any(mask(:))
        metric.r = NaN;
        metric.p_circ = NaN;
        metric.p_entry = NaN;
        metric.CI = [NaN NaN];
        return
    end

    [metric.r, metric.p_circ] = perm_test_circshift(A, B, mask, Nperm);
    [~, metric.p_entry] = perm_test_entry_shuffle(A, B, mask, Nperm);
    [~, metric.CI] = bootstrap_similarity_scalar(A, B, mask, Nboot);
end

function interaction = compute_ssim_region_interaction(dros_struct, marm_struct, comparison_mode, Nperm, peripheral_ssim, central_ssim, ssim_debug)
% Test whether cross-species SSIM is stronger in peripheral than central.
    peripheral_idx = find(strcmp({dros_struct.name}, 'peripheral'), 1);
    central_idx = find(strcmp({dros_struct.name}, 'central'), 1);
    if isempty(peripheral_idx) || isempty(central_idx)
        error('SSIM interaction requires both peripheral and central regions.');
    end

    Xp = dros_struct(peripheral_idx).scaled;
    Yp = marm_struct(peripheral_idx).scaled;
    Xc = dros_struct(central_idx).scaled;
    Yc = marm_struct(central_idx).scaled;

    mask_p = get_matrix_pair_mask(Xp, Yp);
    mask_c = get_matrix_pair_mask(Xc, Yc);

    old_interaction_peripheral_ssim = matrix_ssim(Xp, Yp, mask_p);
    old_interaction_central_ssim = matrix_ssim(Xc, Yc, mask_c);
    old_interaction_diff = old_interaction_peripheral_ssim - old_interaction_central_ssim;

    if nargin >= 7 && ~isempty(ssim_debug)
        print_ssim_path_debug('peripheral', Xp, Yp, mask_p, old_interaction_peripheral_ssim, ssim_debug);
        print_ssim_path_debug('central', Xc, Yc, mask_c, old_interaction_central_ssim, ssim_debug);
    end

    % Use the already reported per-region SSIM values for the observed
    % interaction. This guarantees that PeripheralMinusCentralSSIM is exactly
    % the difference between the values printed in the main comparison table.
    if nargin < 5 || isempty(peripheral_ssim)
        peripheral_ssim = matrix_ssim(Xp, Yp, mask_p);
    end
    if nargin < 6 || isempty(central_ssim)
        central_ssim = matrix_ssim(Xc, Yc, mask_c);
    end
    observed_diff = peripheral_ssim - central_ssim;

    circ_null = nan(Nperm, 1);
    entry_null = nan(Nperm, 1);
    for perm_ind = 1:Nperm
        Xp_circ = circular_shift_matrix_randomly(Xp);
        Xc_circ = circular_shift_matrix_randomly(Xc);
        circ_null(perm_ind) = matrix_ssim(Xp_circ, Yp, mask_p) - matrix_ssim(Xc_circ, Yc, mask_c);

        Xp_entry = shuffle_matrix_entries_within_mask(Xp, mask_p);
        Xc_entry = shuffle_matrix_entries_within_mask(Xc, mask_c);
        entry_null(perm_ind) = matrix_ssim(Xp_entry, Yp, mask_p) - matrix_ssim(Xc_entry, Yc, mask_c);
    end

    p_circ_greater = (sum(circ_null >= observed_diff) + 1) / (Nperm + 1);
    p_entry_greater = (sum(entry_null >= observed_diff) + 1) / (Nperm + 1);
    p_circ_two_sided = (sum(abs(circ_null) >= abs(observed_diff)) + 1) / (Nperm + 1);
    p_entry_two_sided = (sum(abs(entry_null) >= abs(observed_diff)) + 1) / (Nperm + 1);

    interaction.table = table( ...
        string(comparison_mode), peripheral_ssim, central_ssim, observed_diff, ...
        old_interaction_peripheral_ssim, old_interaction_central_ssim, old_interaction_diff, ...
        p_circ_greater, p_entry_greater, p_circ_two_sided, p_entry_two_sided, Nperm, ...
        'VariableNames', {'ComparisonMode','PeripheralSSIM','CentralSSIM','PeripheralMinusCentralSSIM', ...
        'OldInteractionPeripheralSSIM','OldInteractionCentralSSIM','OldInteractionPeripheralMinusCentralSSIM', ...
        'PCircGreater','PEntryGreater','PCircTwoSided','PEntryTwoSided','NPerm'});
    interaction.mode = comparison_mode;
    interaction.observed_diff = observed_diff;
    interaction.peripheral_ssim = peripheral_ssim;
    interaction.central_ssim = central_ssim;
    interaction.old_interaction_peripheral_ssim = old_interaction_peripheral_ssim;
    interaction.old_interaction_central_ssim = old_interaction_central_ssim;
    interaction.old_interaction_diff = old_interaction_diff;
    interaction.circ_null = circ_null;
    interaction.entry_null = entry_null;

    fprintf('  %s: SSIM peripheral=%.3f, central=%.3f, peripheral-central=%.3f, p_circ=%.4f, p_entry=%.4f\n', ...
        comparison_mode, peripheral_ssim, central_ssim, observed_diff, p_circ_greater, p_entry_greater);
    fprintf('    Old interaction recomputation: peripheral=%.12f, central=%.12f, peripheral-central=%.12f\n', ...
        old_interaction_peripheral_ssim, old_interaction_central_ssim, old_interaction_diff);
end

function mask = get_matrix_pair_mask(A, B)
    mask = ~isnan(A) & ~isnan(B);
    if all(~mask(:))
        mask = isfinite(A) & isfinite(B);
    end
end

function print_ssim_path_debug(region_name, X_interaction, Y_interaction, mask_interaction, ssim_interaction, ssim_debug)
    debug_idx = find(strcmp({ssim_debug.region}, region_name), 1);
    if isempty(debug_idx) || isempty(ssim_debug(debug_idx).X)
        fprintf('    SSIM path debug [%s]: no main-loop debug matrix stored.\n', region_name);
        return
    end

    X_main = ssim_debug(debug_idx).X;
    Y_main = ssim_debug(debug_idx).Y;
    mask_main = ssim_debug(debug_idx).mask;
    ssim_main = ssim_debug(debug_idx).ssim;

    if ~isequal(size(X_main), size(X_interaction)) || ~isequal(size(Y_main), size(Y_interaction))
        fprintf('    SSIM path debug [%s]: size differs. main X=%dx%d, interaction X=%dx%d; main Y=%dx%d, interaction Y=%dx%d\n', ...
            region_name, size(X_main,1), size(X_main,2), size(X_interaction,1), size(X_interaction,2), ...
            size(Y_main,1), size(Y_main,2), size(Y_interaction,1), size(Y_interaction,2));
        return
    end

    max_x_diff = max(abs(X_main(:) - X_interaction(:)), [], 'omitnan');
    max_y_diff = max(abs(Y_main(:) - Y_interaction(:)), [], 'omitnan');
    mask_diff = nnz(mask_main ~= mask_interaction);
    ssim_main_recomputed = matrix_ssim(X_main, Y_main, mask_main);

    fprintf('    SSIM path debug [%s]: main stored=%.12f, main recomputed=%.12f, interaction recomputed=%.12f, max|dX|=%.3g, max|dY|=%.3g, mask diff=%d\n', ...
        region_name, ssim_main, ssim_main_recomputed, ssim_interaction, max_x_diff, max_y_diff, mask_diff);
end

function q = bh_fdr(p)
% Benjamini-Hochberg FDR correction, preserving NaNs and original order.
    q = nan(size(p));
    finite_mask = isfinite(p);
    p_finite = p(finite_mask);
    m = numel(p_finite);
    if m == 0
        return
    end

    [p_sorted, sort_idx] = sort(p_finite(:), 'ascend');
    q_sorted = p_sorted .* m ./ (1:m)';
    q_sorted = flipud(cummin(flipud(q_sorted)));
    q_sorted(q_sorted > 1) = 1;

    q_finite = nan(m, 1);
    q_finite(sort_idx) = q_sorted;
    q(finite_mask) = reshape(q_finite, size(p_finite));
end

function [miCurve, miSem] = build_region_mi_curve(miMatrix, sourceRows, targetLength)
% Average pair-level MI curves and resize to the normalized CoI matrix length.
    miCurve = nan(targetLength, 1);
    miSem = nan(targetLength, 1);
    if isempty(miMatrix)
        return
    end

    validRows = sourceRows(sourceRows >= 1 & sourceRows <= size(miMatrix, 1));
    if isempty(validRows)
        return
    end

    regionMI = mean(miMatrix(validRows, :), 2, 'omitnan');
    regionSem = std(miMatrix(validRows, :), 0, 2, 'omitnan') ./ sqrt(sum(isfinite(miMatrix(validRows, :)), 2));
    miCurve = interpolate_vector_to_length(regionMI, targetLength);
    miSem = interpolate_vector_to_length(regionSem, targetLength);
end

function yOut = interpolate_vector_to_length(yIn, targetLength)
    yIn = yIn(:);
    yOut = nan(targetLength, 1);
    valid = isfinite(yIn);
    if nnz(valid) < 2
        return
    end

    sourceX = linspace(1, targetLength, numel(yIn));
    targetX = 1:targetLength;
    yOut = interp1(sourceX(valid), yIn(valid), targetX, 'linear', 'extrap')';
end

function marmMI = load_marmoset_mi_curves(basefold, target_n, target_n_all, target_pre, marmoset_stim_onset_ms)
% Load Figure 2e marmoset MI means prepared by plot_marmoset_mi_figure2e.
    csvFile = fullfile(basefold, 'Results', 'Marmoset', 'Marmoset_MI', 'marmoset_figure2e_mi.csv');
    regionsLocal = {'peripheral', 'central'};
    for idx = 1:numel(regionsLocal)
        marmMI.(regionsLocal{idx}).post = nan(target_n, 1);
        marmMI.(regionsLocal{idx}).post_sem = nan(target_n, 1);
        marmMI.(regionsLocal{idx}).with_prestim = nan(target_n_all, 1);
        marmMI.(regionsLocal{idx}).with_prestim_sem = nan(target_n_all, 1);
    end
    marmMI.time_post = (0:target_n-1)';
    onsetIdx = round(marmoset_stim_onset_ms);
    preStart = max(1, onsetIdx - target_pre);
    marmMI.time_with_prestim = ((preStart:(preStart + target_n_all - 1))' - onsetIdx);

    if ~exist(csvFile, 'file')
        warning('Marmoset MI CSV not found: %s. MI marginals for marmoset will be NaN.', csvFile);
        return
    end

    T = readtable(csvFile);
    peripheralMI = get_table_vector(T, 'TemporalMean_Peripheral');
    peripheralSEM = compute_table_sem(T, {'Temporal_Kr', 'Temporal_Go', 'Temporal_Fr'});
    centralMI = get_table_vector(T, 'FrontalMean_Central');
    centralSEM = compute_table_sem(T, {'Frontal_Kr', 'Frontal_Go', 'Frontal_Fr'});

    postStart = onsetIdx;
    postIdx = postStart:(postStart + target_n - 1);
    preIdx = preStart:(preStart + target_n_all - 1);

    marmMI.peripheral.post = take_or_interpolate_mi_window(peripheralMI, postIdx, target_n);
    marmMI.peripheral.post_sem = take_or_interpolate_mi_window(peripheralSEM, postIdx, target_n);
    marmMI.central.post = take_or_interpolate_mi_window(centralMI, postIdx, target_n);
    marmMI.central.post_sem = take_or_interpolate_mi_window(centralSEM, postIdx, target_n);
    marmMI.peripheral.with_prestim = take_or_interpolate_mi_window(peripheralMI, preIdx, target_n_all);
    marmMI.peripheral.with_prestim_sem = take_or_interpolate_mi_window(peripheralSEM, preIdx, target_n_all);
    marmMI.central.with_prestim = take_or_interpolate_mi_window(centralMI, preIdx, target_n_all);
    marmMI.central.with_prestim_sem = take_or_interpolate_mi_window(centralSEM, preIdx, target_n_all);
end

function y = get_table_vector(T, varName)
    if ~ismember(varName, T.Properties.VariableNames)
        error('Could not find column %s in marmoset MI CSV.', varName);
    end
    y = T.(varName);
    y = y(:);
end

function semCurve = compute_table_sem(T, varNames)
    values = nan(height(T), numel(varNames));
    for colIdx = 1:numel(varNames)
        if ~ismember(varNames{colIdx}, T.Properties.VariableNames)
            error('Could not find column %s in marmoset MI CSV.', varNames{colIdx});
        end
        values(:, colIdx) = T.(varNames{colIdx});
    end
    semCurve = std(values, 0, 2, 'omitnan') ./ sqrt(sum(isfinite(values), 2));
end

function yWindow = take_or_interpolate_mi_window(y, idx, targetLength)
    yWindow = nan(targetLength, 1);
    validIdx = idx(idx >= 1 & idx <= numel(y));
    if numel(validIdx) == targetLength
        yWindow = y(validIdx);
    elseif ~isempty(validIdx)
        yWindow = interpolate_vector_to_length(y(validIdx), targetLength);
    end
end

function fig = plot_cross_species_matrices_with_mi(drosStruct, marmStruct, regionsLocal, labels, drosTime, marmTime, titleText)
    marmosetLabels = map_marmoset_region_labels(regionsLocal);
    fig = figure('Color', 'w', 'Position', [60 40 1600 1150], 'Renderer', 'painters');
    sgtitle(titleText, 'Interpreter', 'none');

    % Manual axes positions keep the four CoI matrices exactly the same size.
    % Each MI axis is then positioned from its paired matrix axis, so it uses
    % the same time vector and cannot drift out of alignment with the CoI.
    matrixSide = 0.255;
    matrixX = [0.235, 0.565];
    matrixY = [0.565, 0.170]; % top row, bottom row
    gap = 0.014;
    miWidth = 0.075;
    miHeight = 0.075;

    panelSpecs = struct( ...
        'matrixPos', { ...
            [matrixX(1), matrixY(1), matrixSide, matrixSide], ...
            [matrixX(2), matrixY(1), matrixSide, matrixSide], ...
            [matrixX(1), matrixY(2), matrixSide, matrixSide], ...
            [matrixX(2), matrixY(2), matrixSide, matrixSide]}, ...
        'hSide', {'top', 'top', 'bottom', 'bottom'}, ...
        'vSide', {'left', 'right', 'left', 'right'}, ...
        'matrixData', {drosStruct(1).scaled, drosStruct(2).scaled, marmStruct(1).scaled, marmStruct(2).scaled}, ...
        'miCurve', {drosStruct(1).mi, drosStruct(2).mi, marmStruct(1).mi, marmStruct(2).mi}, ...
        'miSem', {drosStruct(1).mi_sem, drosStruct(2).mi_sem, marmStruct(1).mi_sem, marmStruct(2).mi_sem}, ...
        'timeAxis', {drosTime, drosTime, marmTime, marmTime}, ...
        'panelTitle', { ...
            sprintf('Drosophila - %s', char(labels(1))), ...
            sprintf('Drosophila - %s', char(labels(2))), ...
            sprintf('Marmoset - %s', marmosetLabels{1}), ...
            sprintf('Marmoset - %s', marmosetLabels{2})}, ...
        'miColor', {[0.10 0.35 0.85], [0.10 0.35 0.85], [0.85 0.25 0.10], [0.85 0.25 0.10]});

    drosMiLimits = get_common_mi_limits(panelSpecs(1:2));
    marmMiLimits = get_common_mi_limits(panelSpecs(3:4));
    for panelInd = 1:numel(panelSpecs)
        if panelInd <= 2
            miLimits = drosMiLimits;
        else
            miLimits = marmMiLimits;
        end
        plot_outer_mi_matrix_panel(fig, panelSpecs(panelInd), gap, miWidth, miHeight, miLimits);
    end

    axCb = axes('Parent', fig, 'Visible', 'off', 'Position', [0.935 0.170 0.001 0.650]);
    colormap(axCb, redblue(256));
    caxis(axCb, [-2 2]);
    cb = colorbar(axCb);
    cb.Position = [0.940 0.315 0.009 0.385];
    cb.Ticks = [-2 0 2];
    cb.TickLabels = {'-2 Synergistic', '0', '2 Redundant'};
    cb.Label.String = 'Normalized CoI (bits)';
    cb.Label.Rotation = 270;
    cb.Label.VerticalAlignment = 'bottom';
    cb.FontSize = 12;
    cb.Label.FontSize = 14;
    cb.LineWidth = 0.45;
end

function plot_outer_mi_matrix_panel(fig, panelSpec, gap, miWidth, miHeight, miLimits)
    matrixPos = panelSpec.matrixPos;

    [matrixData, miCurve, miSem, timeAxis] = prepare_matrix_and_mi( ...
        panelSpec.matrixData, panelSpec.miCurve, panelSpec.miSem, panelSpec.timeAxis);

    axMat = axes('Parent', fig, 'Position', matrixPos);
    lock_axes_position(axMat);
    imagesc(axMat, timeAxis, timeAxis, matrixData);
    set(axMat, 'YDir', 'normal', 'FontSize', 9, 'TickDir', 'out', ...
        'Layer', 'top', 'LineWidth', 0.45, 'Box', 'on');
    axis(axMat, 'square');
    xlim(axMat, [min(timeAxis), max(timeAxis)]);
    ylim(axMat, [min(timeAxis), max(timeAxis)]);
    colormap(axMat, redblue(256));
    caxis(axMat, [-2 2]);
    xlabel(axMat, 'Time (ms)');
    ylabel(axMat, 'Time (ms)');
    title(axMat, panelSpec.panelTitle, 'Interpreter', 'none');

    drawnow;
    matrixPlotBox = get_square_plot_box_position(axMat);
    matrixOuterBox = get_axes_outer_position_normalized(axMat);
    hPos = get_horizontal_mi_position(matrixPlotBox, matrixOuterBox, panelSpec.hSide, gap, miHeight);
    vPos = get_vertical_mi_position(matrixPlotBox, panelSpec.vSide, gap, miWidth);

    axH = axes('Parent', fig, 'Position', hPos);
    lock_axes_position(axH);
    plot_horizontal_mi(axH, timeAxis, miCurve, miSem, panelSpec.miColor, panelSpec.hSide, miLimits);

    axV = axes('Parent', fig, 'Position', vPos);
    lock_axes_position(axV);
    plot_vertical_mi(axV, timeAxis, miCurve, miSem, panelSpec.miColor, panelSpec.vSide, miLimits);

    check_marginal_alignment(matrixPlotBox, get(axH, 'Position'), get(axV, 'Position'), panelSpec.panelTitle);
end

function miLimits = get_common_mi_limits(panelSpecs)
    allValues = [];
    for panelInd = 1:numel(panelSpecs)
        miCurve = panelSpecs(panelInd).miCurve(:);
        miSem = panelSpecs(panelInd).miSem(:);
        if isempty(miSem)
            miSem = zeros(size(miCurve));
        end
        if numel(miSem) ~= numel(miCurve)
            miSem = interpolate_vector_to_length(miSem, numel(miCurve));
        end
        allValues = [allValues; miCurve - miSem; miCurve + miSem]; %#ok<AGROW>
    end

    allValues = allValues(isfinite(allValues));
    if isempty(allValues)
        miLimits = [0 1];
        return
    end

    lo = min(allValues);
    hi = max(allValues);
    if lo == hi
        pad = max(abs(lo) * 0.1, 1);
    else
        pad = (hi - lo) * 0.08;
    end
    miLimits = [lo - pad, hi + pad];
end

function [matrixData, miCurve, miSem, timeAxis] = prepare_matrix_and_mi(matrixData, miCurve, miSem, timeAxis)
    n = size(matrixData, 1);
    if numel(timeAxis) ~= n
        timeAxis = linspace(min(timeAxis), max(timeAxis), n);
    end
    timeAxis = timeAxis(:)';
    miCurve = interpolate_vector_to_length(miCurve, n);
    miSem = interpolate_vector_to_length(miSem, n);
end

function pos = get_horizontal_mi_position(matrixPos, matrixOuterBox, sideName, gap, miHeight)
    % Objective alignment: matrixPos is the actual square CoI image plot box,
    % not the wider axes box that MATLAB uses before applying axis square.
    if strcmp(sideName, 'top')
        pos = [matrixPos(1), matrixPos(2) + matrixPos(4) + gap, matrixPos(3), miHeight];
    else
        % For bottom marginals, place the MI below the full matrix annotation
        % extent so x tick labels and "Time (ms)" do not collide with the MI.
        pos = [matrixPos(1), matrixOuterBox(2) - gap - miHeight, matrixPos(3), miHeight];
    end
end

function pos = get_vertical_mi_position(matrixPos, sideName, gap, miWidth)
    if strcmp(sideName, 'left')
        leftGap = gap + 0.012;
        pos = [matrixPos(1) - leftGap - miWidth, matrixPos(2), miWidth, matrixPos(4)];
    else
        pos = [matrixPos(1) + matrixPos(3) + gap, matrixPos(2), miWidth, matrixPos(4)];
    end
end

function plot_horizontal_mi(ax, timeAxis, miCurve, miSem, miColor, sideName, miLimits)
    miLower = miCurve - miSem;
    miUpper = miCurve + miSem;

    hold(ax, 'on');
    plot_sem_patch(ax, timeAxis, miLower, miUpper, miColor, false);
    plot(ax, timeAxis, miCurve, 'Color', miColor, 'LineWidth', 0.75);
    xlim(ax, [min(timeAxis), max(timeAxis)]);
    ylim(ax, miLimits);
    if strcmp(sideName, 'bottom')
        set(ax, 'YDir', 'reverse', 'XAxisLocation', 'bottom');
    else
        set(ax, 'XAxisLocation', 'top');
    end
    ylabel(ax, 'MI (bits)');
    set(ax, 'XTick', [], 'Box', 'on', 'TickDir', 'out', ...
        'FontSize', 8, 'LineWidth', 0.45, ...
        'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15]);
    hold(ax, 'off');
end

function plot_vertical_mi(ax, timeAxis, miCurve, miSem, miColor, sideName, miLimits)
    miLower = miCurve - miSem;
    miUpper = miCurve + miSem;

    hold(ax, 'on');
    plot_sem_patch(ax, timeAxis, miLower, miUpper, miColor, true);
    plot(ax, miCurve, timeAxis, 'Color', miColor, 'LineWidth', 0.75);
    ylim(ax, [min(timeAxis), max(timeAxis)]);
    xlim(ax, miLimits);
    if strcmp(sideName, 'left')
        set(ax, 'XDir', 'reverse');
    else
        set(ax, 'YAxisLocation', 'right');
    end
    xlabel(ax, 'MI (bits)');
    set(ax, 'YTickLabel', [], 'Box', 'on', 'TickDir', 'out', ...
        'FontSize', 8, 'LineWidth', 0.45, ...
        'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15]);
    hold(ax, 'off');
end

function plotBox = get_square_plot_box_position(ax)
    fig = ancestor(ax, 'figure');
    oldAxUnits = get(ax, 'Units');
    oldFigUnits = get(fig, 'Units');
    set(ax, 'Units', 'pixels');
    set(fig, 'Units', 'pixels');

    axPix = get(ax, 'Position');
    figPix = get(fig, 'Position');
    sidePix = min(axPix(3), axPix(4));
    plotPix = [ ...
        axPix(1) + (axPix(3) - sidePix) / 2, ...
        axPix(2) + (axPix(4) - sidePix) / 2, ...
        sidePix, ...
        sidePix];

    plotBox = [ ...
        plotPix(1) / figPix(3), ...
        plotPix(2) / figPix(4), ...
        plotPix(3) / figPix(3), ...
        plotPix(4) / figPix(4)];

    set(ax, 'Units', oldAxUnits);
    set(fig, 'Units', oldFigUnits);
end

function outerBox = get_axes_outer_position_normalized(ax)
    oldAxUnits = get(ax, 'Units');
    set(ax, 'Units', 'normalized');
    outerBox = get(ax, 'OuterPosition');
    set(ax, 'Units', oldAxUnits);
end

function lock_axes_position(ax)
    if isprop(ax, 'PositionConstraint')
        set(ax, 'PositionConstraint', 'innerposition');
    elseif isprop(ax, 'ActivePositionProperty')
        set(ax, 'ActivePositionProperty', 'position');
    end
end

function check_marginal_alignment(matrixPos, hPos, vPos, panelTitle)
    tol = 1e-10;
    hAligned = abs(hPos(1) - matrixPos(1)) < tol && abs(hPos(3) - matrixPos(3)) < tol;
    vAligned = abs(vPos(2) - matrixPos(2)) < tol && abs(vPos(4) - matrixPos(4)) < tol;
    if ~(hAligned && vAligned)
        warning('MI marginal axes are not aligned for %s.', panelTitle);
    end
end

function labelsOut = map_marmoset_region_labels(regionsLocal)
    labelsOut = cell(size(regionsLocal));
    for idx = 1:numel(regionsLocal)
        switch char(regionsLocal{idx})
            case 'peripheral'
                labelsOut{idx} = 'Temporal';
            case 'central'
                labelsOut{idx} = 'Frontal';
            otherwise
                labelsOut{idx} = strrep(char(regionsLocal{idx}), '_', '-');
        end
    end
end

function plot_sem_patch(ax, timeAxis, lowerCurve, upperCurve, colorValue, sideways)
    valid = isfinite(timeAxis(:)) & isfinite(lowerCurve(:)) & isfinite(upperCurve(:));
    if nnz(valid) < 2
        return
    end
    t = timeAxis(valid);
    lo = lowerCurve(valid);
    hi = upperCurve(valid);

    if sideways
        patch(ax, [lo(:); flipud(hi(:))], [t(:); flipud(t(:))], colorValue, ...
            'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    else
        patch(ax, [t(:); flipud(t(:))], [lo(:); flipud(hi(:))], colorValue, ...
            'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

function shifted = circular_shift_matrix_randomly(A)
    n_rows = size(A, 1);
    n_cols = size(A, 2);
    shifted = circshift(A, randi(n_rows), 1);
    shifted = circshift(shifted, randi(n_cols), 2);
end

function shuffled = shuffle_matrix_entries_within_mask(A, mask)
    shuffled = A;
    vals = A(mask);
    valid = isfinite(vals);
    finite_vals = vals(valid);
    if numel(finite_vals) < 2
        return
    end
    vals(valid) = finite_vals(randperm(numel(finite_vals)));
    shuffled(mask) = vals;
end

function plot_ssim_interaction_nulls(interaction_post, interaction_with_prestim, out_dir, activity_tag, condition)
    interactions = {interaction_post, interaction_with_prestim};
    fig = figure('Color', 'w', 'Position', [100 100 1050 440], 'Renderer', 'painters');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for interaction_ind = 1:numel(interactions)
        interaction = interactions{interaction_ind};
        nexttile;
        hold on
        histogram(interaction.circ_null, 40, 'FaceColor', [0.25 0.45 0.75], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
        histogram(interaction.entry_null, 40, 'FaceColor', [0.85 0.35 0.15], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        xline(interaction.observed_diff, 'k-', 'LineWidth', 0.75);
        xline(0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.5);
        xlabel('Peripheral SSIM - Central SSIM');
        ylabel('Permutation count');
        title(strrep(interaction.mode, '_', '\_'));
        legend({'Circular shift null', 'Entry shuffle null', 'Observed'}, 'Box', 'off', 'Location', 'best');
        grid on
        box off
        set(gca, 'LineWidth', 0.45);
        hold off
    end

    sgtitle('SSIM Region Interaction: Cross-Species Similarity Difference');
    saveas(fig, fullfile(out_dir, sprintf('ssim_region_interaction_%s_%s.png', activity_tag, condition)));
    saveas(fig, fullfile(out_dir, sprintf('ssim_region_interaction_%s_%s.fig', activity_tag, condition)));
    print_editable_svg(fig, fullfile(out_dir, sprintf('ssim_region_interaction_%s_%s.svg', activity_tag, condition)));
    close(fig);
end

function plot_ssim_region_interaction_bar(results_table, interaction_results, comparison_mode)
    peripheral_idx = find(strcmp(results_table.Region, 'peripheral'), 1);
    central_idx = find(strcmp(results_table.Region, 'central'), 1);
    interaction_idx = find(interaction_results.ComparisonMode == comparison_mode, 1);

    values = [ ...
        results_table.ssim(peripheral_idx), ...
        results_table.ssim(central_idx), ...
        interaction_results.PeripheralMinusCentralSSIM(interaction_idx)];
    labels = {'Peripheral-temporal', 'Central-frontal', 'Interaction'};
    colors = [0.30 0.55 0.80; 0.80 0.45 0.25; 0.35 0.35 0.35];

    b = bar(values, 'FaceColor', 'flat', 'LineWidth', 0.45);
    b.CData = colors;
    hold on
    yline(0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.45);

    xticks(1:3);
    xticklabels(labels);
    xtickangle(20);
    ylabel('SSIM / SSIM difference');
    title('SSIM Region Summary');
    grid on
    box off
    set(gca, 'LineWidth', 0.45);

    p_values = [ ...
        results_table.ssim_q_circ(peripheral_idx), ...
        results_table.ssim_q_circ(central_idx), ...
        interaction_p_value_for_plot(interaction_results, interaction_idx)];

    y_min = min([0, values], [], 'omitnan');
    y_max = max([0, values], [], 'omitnan');
    y_range = max(y_max - y_min, eps);
    ylim([y_min - 0.12*y_range, y_max + 0.22*y_range]);

    for value_idx = 1:3
        star_txt = p_to_star_local(p_values(value_idx));
        if isempty(star_txt)
            continue
        end
        text(value_idx, values(value_idx) + 0.05*y_range, star_txt, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 16, ...
            'FontWeight', 'bold');
    end
    hold off
end

function p_val = interaction_p_value_for_plot(interaction_results, interaction_idx)
    if ismember('QCircGreater', interaction_results.Properties.VariableNames) && ...
            isfinite(interaction_results.QCircGreater(interaction_idx))
        p_val = interaction_results.QCircGreater(interaction_idx);
    else
        p_val = interaction_results.PCircGreater(interaction_idx);
    end
end

function star_txt = p_to_star_local(p_val)
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

% -----------------------------------------------------------------
% STEP 5: Save results to CSV
% -----------------------------------------------------------------
outcsv = fullfile(out_dir, sprintf('drosophila_marmoset_comparison_%s_%s.csv', activity_tag, condition));
writetable(results, outcsv);
fprintf('Results saved to: %s\n', outcsv);

outcsv_pre = fullfile(out_dir, sprintf('drosophila_marmoset_comparison_with_prestim_%s_%s.csv', activity_tag, condition));
writetable(results_with_prestim, outcsv_pre);
fprintf('Results (with pre-stim) saved to: %s\n', outcsv_pre);

outcsv_interaction = fullfile(out_dir, sprintf('drosophila_marmoset_ssim_region_interaction_%s_%s.csv', activity_tag, condition));
writetable(interaction_results, outcsv_interaction);
fprintf('SSIM region interaction results saved to: %s\n', outcsv_interaction);

% -----------------------------------------------------------------
% STEP 6: Generate figures
% -----------------------------------------------------------------
fprintf('Generating comparison figures...\n');

% Figure 1: Observed Pearson r, SSIM per region, and SSIM interaction
fig1 = figure('Position', [100 100 1350 400], 'Renderer', 'painters');
subplot(1, 3, 1);
bar(results.r_pearson, 'FaceColor', [0.2 0.5 0.8]);
hold on;
errorbar(1:n_regions, results.r_pearson, results.r_pearson - results.CI_low, ...
    results.CI_high - results.r_pearson, 'k.', 'MarkerSize', 1, 'CapSize', 5, 'LineWidth', 0.55);
xticks(1:n_regions);
xticklabels(results.Label);
ylabel('Pearson r');
title('Drosophila-Marmoset CoI Correlation');
grid on; ylim([-1 1]); set(gca, 'LineWidth', 0.45);
sig_circ = results.p_circ < 0.05;
for i = 1:n_regions
    if sig_circ(i)
        text(i, results.r_pearson(i) + 0.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 16);
    end
end
hold off;

subplot(1, 3, 2);
bar(results.ssim, 'FaceColor', [0.8 0.5 0.2]);
hold on;
errorbar(1:n_regions, results.ssim, results.ssim - results.ssim_CI_low, ...
    results.ssim_CI_high - results.ssim, 'k.', 'MarkerSize', 1, 'CapSize', 5, 'LineWidth', 0.55);
xticks(1:n_regions);
xticklabels(results.Label);
ylabel('SSIM');
title('Structural Similarity Index');
grid on; ylim([0 1]); set(gca, 'LineWidth', 0.45);
sig_circ_ssim = results.ssim_q_circ < 0.05;
for i = 1:n_regions
    if sig_circ_ssim(i)
        text(i, results.ssim(i) + 0.08, '*', 'HorizontalAlignment', 'center', 'FontSize', 16);
    end
end
hold off;

subplot(1, 3, 3);
plot_ssim_region_interaction_bar(results, interaction_results, "post_only");

saveas(fig1, fullfile(out_dir, sprintf('cross_species_metrics_%s_%s.png', activity_tag, condition)));
print_editable_svg(fig1, fullfile(out_dir, sprintf('cross_species_metrics_%s_%s.svg', activity_tag, condition)));
close(fig1);

% Figure 1b: Observed Pearson r, SSIM per region, and SSIM interaction (with pre-stim)
fig1_pre = figure('Position', [100 100 1350 400], 'Renderer', 'painters');
subplot(1, 3, 1);
bar(results_with_prestim.r_pearson, 'FaceColor', [0.2 0.5 0.8]);
hold on;
errorbar(1:n_regions, results_with_prestim.r_pearson, results_with_prestim.r_pearson - results_with_prestim.CI_low, ...
    results_with_prestim.CI_high - results_with_prestim.r_pearson, 'k.', 'MarkerSize', 1, 'CapSize', 5, 'LineWidth', 0.55);
xticks(1:n_regions);
xticklabels(results_with_prestim.Label);
ylabel('Pearson r');
title('Drosophila-Marmoset CoI Correlation (With Pre-Stim)');
grid on; ylim([-1 1]); set(gca, 'LineWidth', 0.45);
sig_circ_pre = results_with_prestim.p_circ < 0.05;
for i = 1:n_regions
    if sig_circ_pre(i)
        text(i, results_with_prestim.r_pearson(i) + 0.15, '*', 'HorizontalAlignment', 'center', 'FontSize', 16);
    end
end
hold off;

subplot(1, 3, 2);
bar(results_with_prestim.ssim, 'FaceColor', [0.8 0.5 0.2]);
hold on;
errorbar(1:n_regions, results_with_prestim.ssim, results_with_prestim.ssim - results_with_prestim.ssim_CI_low, ...
    results_with_prestim.ssim_CI_high - results_with_prestim.ssim, 'k.', 'MarkerSize', 1, 'CapSize', 5, 'LineWidth', 0.55);
xticks(1:n_regions);
xticklabels(results_with_prestim.Label);
ylabel('SSIM');
title('Structural Similarity Index (With Pre-Stim)');
grid on; ylim([0 1]); set(gca, 'LineWidth', 0.45);
sig_circ_ssim_pre = results_with_prestim.ssim_q_circ < 0.05;
for i = 1:n_regions
    if sig_circ_ssim_pre(i)
        text(i, results_with_prestim.ssim(i) + 0.08, '*', 'HorizontalAlignment', 'center', 'FontSize', 16);
    end
end
hold off;

subplot(1, 3, 3);
plot_ssim_region_interaction_bar(results_with_prestim, interaction_results, "with_prestim");

saveas(fig1_pre, fullfile(out_dir, sprintf('cross_species_metrics_with_prestim_%s_%s.png', activity_tag, condition)));
print_editable_svg(fig1_pre, fullfile(out_dir, sprintf('cross_species_metrics_with_prestim_%s_%s.svg', activity_tag, condition)));
close(fig1_pre);

% Figure 2: Drosophila vs marmoset matrices (scaled) for each region
fig2 = figure('Position', [100 100 1400 600], 'Renderer', 'painters');
sgtitle(sprintf('Scaled CoI Matrices: Drosophila vs Marmoset (%s, %s)', activity_tag, condition));
tiledlayout(2, n_regions, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for r = 1:n_regions
    % Drosophila
    nexttile(r);
    X = dros_struct(r).scaled;
    imagesc(X);
    set(gca, 'YDir', 'normal', 'LineWidth', 0.45);
    axis square; cbTmp = colorbar; cbTmp.LineWidth = 0.45; cbTmp.FontSize = 10; colormap(gca, redblue(256));
    caxis([-2 2]);
    title(sprintf('Drosophila - %s', results.Label(r)));
    
    % Marmoset
    nexttile(n_regions + r);
    Y = marm_struct(r).scaled;
    imagesc(Y);
    set(gca, 'YDir', 'normal', 'LineWidth', 0.45);
    axis square; cbTmp = colorbar; cbTmp.LineWidth = 0.45; cbTmp.FontSize = 10; colormap(gca, redblue(256));
    caxis([-2 2]);
    title(sprintf('Marmoset - %s', results.Label(r)));
end

saveas(fig2, fullfile(out_dir, sprintf('cross_species_matrices_%s_%s.png', activity_tag, condition)));
print_editable_svg(fig2, fullfile(out_dir, sprintf('cross_species_matrices_%s_%s.svg', activity_tag, condition)));
close(fig2);

% Figure 2b: Drosophila vs marmoset matrices (scaled) including pre-stim
fig2_pre = figure('Position', [100 100 1400 600], 'Renderer', 'painters');
sgtitle(sprintf('Scaled CoI Matrices (With Pre-Stim): Drosophila vs Marmoset (%s, %s)', activity_tag, condition));
tiledlayout(2, n_regions, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for r = 1:n_regions
    % Drosophila (with pre-stim)
    nexttile(r);
    Xp = dros_struct_prestim(r).scaled;
    imagesc(Xp);
    set(gca, 'YDir', 'normal', 'LineWidth', 0.45);
    axis square; cbTmp = colorbar; cbTmp.LineWidth = 0.45; cbTmp.FontSize = 10; colormap(gca, redblue(256));
    caxis([-2 2]);
    title(sprintf('Drosophila (With Pre) - %s', results_with_prestim.Label(r)));

    % Marmoset (with pre-stim)
    nexttile(n_regions + r);
    Yp = marm_struct_prestim(r).scaled;
    imagesc(Yp);
    set(gca, 'YDir', 'normal', 'LineWidth', 0.45);
    axis square; cbTmp = colorbar; cbTmp.LineWidth = 0.45; cbTmp.FontSize = 10; colormap(gca, redblue(256));
    caxis([-2 2]);
    title(sprintf('Marmoset (With Pre) - %s', results_with_prestim.Label(r)));
end

saveas(fig2_pre, fullfile(out_dir, sprintf('cross_species_matrices_with_prestim_%s_%s.png', activity_tag, condition)));
print_editable_svg(fig2_pre, fullfile(out_dir, sprintf('cross_species_matrices_with_prestim_%s_%s.svg', activity_tag, condition)));
close(fig2_pre);

% Figure 2c: Normalized matrices with MI marginal panels.
fig2_mi = plot_cross_species_matrices_with_mi( ...
    dros_struct, marm_struct, regions, results.Label, ...
    dros_time_post, marm_time_post, ...
    sprintf('Scaled CoI Matrices with MI Marginals: Drosophila vs Marmoset (%s, %s)', activity_tag, condition));
saveas(fig2_mi, fullfile(out_dir, sprintf('cross_species_matrices_with_mi_%s_%s.png', activity_tag, condition)));
saveas(fig2_mi, fullfile(out_dir, sprintf('cross_species_matrices_with_mi_%s_%s.fig', activity_tag, condition)));
print_editable_svg(fig2_mi, fullfile(out_dir, sprintf('cross_species_matrices_with_mi_%s_%s.svg', activity_tag, condition)));
close(fig2_mi);

fig2_pre_mi = plot_cross_species_matrices_with_mi( ...
    dros_struct_prestim, marm_struct_prestim, regions, results_with_prestim.Label, ...
    dros_time_prestim, marm_time_prestim, ...
    sprintf('Scaled CoI Matrices with MI Marginals (With Pre-Stim): Drosophila vs Marmoset (%s, %s)', activity_tag, condition));
saveas(fig2_pre_mi, fullfile(out_dir, sprintf('cross_species_matrices_with_prestim_with_mi_%s_%s.png', activity_tag, condition)));
saveas(fig2_pre_mi, fullfile(out_dir, sprintf('cross_species_matrices_with_prestim_with_mi_%s_%s.fig', activity_tag, condition)));
print_editable_svg(fig2_pre_mi, fullfile(out_dir, sprintf('cross_species_matrices_with_prestim_with_mi_%s_%s.svg', activity_tag, condition)));
close(fig2_pre_mi);

% Return both comparisons as one table with mode labels.
results.ComparisonMode = repmat("post_only", n_regions, 1);
results_with_prestim.ComparisonMode = repmat("with_prestim", n_regions, 1);
results = [results; results_with_prestim];
results.Properties.UserData.SSIMRegionInteraction = interaction_results;

fprintf('Figures saved to: %s\n', out_dir);
fprintf('Cross-species comparison complete!\n');

end

%% Local helper functions (adapted from compare_region_patterns.m)

function print_editable_svg(fig_handle, svg_file)
% Keep SVGs editable in Illustrator and avoid heavy default vector strokes.
thin_figure_for_illustrator(fig_handle);
set(fig_handle, 'Renderer', 'painters');
print(fig_handle, svg_file, '-dsvg', '-vector');
end

function thin_figure_for_illustrator(fig_handle)
% MATLAB SVGs can look visually fine but import into Illustrator with thick
% strokes. Thin line-like objects immediately before export.
thinLineWidth = 0.55;
thinAxisWidth = 0.45;

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
        set(axesObjects(axInd), 'LineWidth', thinAxisWidth);
    catch
    end
end

colorbars = findall(fig_handle, 'Type', 'ColorBar');
for cbInd = 1:numel(colorbars)
    try
        set(colorbars(cbInd), 'LineWidth', thinAxisWidth, 'FontSize', max(get(colorbars(cbInd), 'FontSize'), 12));
        colorbars(cbInd).Label.FontSize = max(colorbars(cbInd).Label.FontSize, 14);
    catch
    end
end
end

function [r, pperm, r_null] = perm_test_circshift(A, B, mask, Nperm)
if nargin < 4 || isempty(Nperm); Nperm = 1000; end
if nargin < 3 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
r = matrix_corr(A, B, mask);
n = size(A, 1);
r_null = nan(Nperm, 1);
for t = 1:Nperm
    sh = randi(n);
    Ashift = circshift(A, sh, 1);
    Ashift = circshift(Ashift, sh, 2);
    r_null(t) = matrix_corr(Ashift, B, mask);
end
pperm = (sum(r_null >= r) + 1) / (Nperm + 1);
end

function [r, pperm, r_null] = perm_test_entry_shuffle(A, B, mask, Nperm)
if nargin < 4 || isempty(Nperm); Nperm = 1000; end
if nargin < 3 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
vecA = A(mask); vecB = B(mask);
valid = ~(isnan(vecA) | isnan(vecB));
vecA = vecA(valid); vecB = vecB(valid);
r = corr(vecA(:), vecB(:), 'Rows', 'complete');
n = numel(vecA);
r_null = nan(Nperm, 1);
for t = 1:Nperm
    idx = randperm(n);
    r_null(t) = corr(vecA(idx), vecB, 'Rows', 'complete');
end
pperm = (sum(r_null >= r) + 1) / (Nperm + 1);
end

function [r, p] = matrix_corr(A, B, mask)
if nargin < 3 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
vecA = A(mask); vecB = B(mask);
valid = ~(isnan(vecA) | isnan(vecB));
vecA = vecA(valid); vecB = vecB(valid);
if isempty(vecA)
    r = NaN; p = NaN; return;
end
[r, p] = corr(vecA(:), vecB(:), 'Type', 'Pearson', 'Rows', 'complete');
end

function s = matrix_ssim(A, B, mask)
if nargin < 3 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
A2 = A; B2 = B;
A2(~mask) = 0; B2(~mask) = 0;
v = [A2(mask); B2(mask)];
minv = min(v); maxv = max(v);
if maxv > minv
    An = (A2 - minv) / (maxv - minv);
    Bn = (B2 - minv) / (maxv - minv);
else
    An = zeros(size(A2)); Bn = zeros(size(B2));
end
if exist('ssim', 'file') == 2
    s = ssim(An, Bn);
else
    s = matrix_corr(An, Bn, mask);
end
end

function [p, obs, r_null] = perm_test_circshift_metric(metricFcn, A, B, mask, Nperm)
if nargin < 5 || isempty(Nperm); Nperm = 1000; end
if nargin < 4 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
obs = metricFcn(A, B, mask);
n = size(A, 1);
r_null = nan(Nperm, 1);
for t = 1:Nperm
    sh = randi(n);
    Ashift = circshift(A, sh, 1);
    Ashift = circshift(Ashift, sh, 2);
    r_null(t) = metricFcn(Ashift, B, mask);
end
p = (sum(r_null >= obs) + 1) / (Nperm + 1);
end

function [p, obs, r_null] = perm_test_entry_shuffle_metric(metricFcn, A, B, mask, Nperm)
if nargin < 5 || isempty(Nperm); Nperm = 1000; end
if nargin < 4 || isempty(mask); mask = ~isnan(A) & ~isnan(B); end
obs = metricFcn(A, B, mask);
vecA = A(mask); vecB = B(mask);
valid = ~(isnan(vecA) | isnan(vecB));
vecA = vecA(valid);
n = numel(vecA);
r_null = nan(Nperm, 1);
for t = 1:Nperm
    idx = randperm(n);
    Aprime = A;
    tmp = A(mask);
    tmp(valid) = vecA(idx);
    Aprime(mask) = tmp;
    r_null(t) = metricFcn(Aprime, B, mask);
end
p = (sum(r_null >= obs) + 1) / (Nperm + 1);
end

function [r_mean, CI] = bootstrap_similarity_scalar(A, B, mask, Nboot)
% For comparing two single matrices (not 3D stacks)
if nargin < 4 || isempty(Nboot); Nboot = 1000; end
if isempty(A) || isempty(B)
    r_mean = NaN; CI = [NaN NaN]; return;
end
% Bootstrap by adding random noise proportional to standard deviation
r_boot = nan(Nboot, 1);
for b = 1:Nboot
    % Resample with replacement
    A_noisy = A + randn(size(A)) * nanstd(A(:)) * 0.01;
    B_noisy = B + randn(size(B)) * nanstd(B(:)) * 0.01;
    r_boot(b) = matrix_corr(A_noisy, B_noisy, mask);
end
r_mean = mean(r_boot);
CI = prctile(r_boot, [2.5 97.5]);
end

function [ssim_mean, CI] = bootstrap_ssim_scalar(A, B, mask, Nboot)
% Bootstrap CI for SSIM (two single matrices)
if nargin < 4 || isempty(Nboot); Nboot = 1000; end
if isempty(A) || isempty(B)
    ssim_mean = NaN; CI = [NaN NaN]; return;
end
ssim_boot = nan(Nboot, 1);
for b = 1:Nboot
    A_noisy = A + randn(size(A)) * nanstd(A(:)) * 0.01;
    B_noisy = B + randn(size(B)) * nanstd(B(:)) * 0.01;
    ssim_boot(b) = matrix_ssim(A_noisy, B_noisy, mask);
end
ssim_mean = mean(ssim_boot);
CI = prctile(ssim_boot, [2.5 97.5]);
end

