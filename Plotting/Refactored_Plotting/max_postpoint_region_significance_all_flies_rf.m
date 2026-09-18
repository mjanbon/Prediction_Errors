%% File: max_postpoint_region_significance_all_flies_rf.m
% Group-level analysis of post-point CoI significance across flies.
%
% For each anatomical region, this script counts how many electrode-pair
% instances show significance in more than 50% of the cells inside two
% post-point windows:
%   1) time1 >= -15 ms and time2 >= 0 ms
%   2) time1 >= 0 ms and time2 >= -15 ms
%
% The result is reported as x out of y pair instances per region, pooled
% across all flies.

function stats = max_postpoint_region_significance_all_flies_rf(combine_windows)
    set_default_plotting();

    %% Parameters
    USING_HPC = 0;
    get_elec = 0;
    cutoff = 10;
    timing = -25:1:74;
    % Fraction of cells inside the window required to count a pair as
    % 'significant' (e.g. 0.5 means >50% of cells must be significant)
    frac_sig_threshold = 0.75;
    % When called with an input, that value overrides the default. Default
    % is true (count the union of defined windows to avoid double-counting).
    if nargin < 1 || isempty(combine_windows)
        combine_windows = true;
    end

    %% Load participant info
    [basefold, datatype, all_con, ~, ~, participants, ~, ~, activity_tag, ~, ~, ~, ~, ~, ~, ~, ~] = Max_get_param(USING_HPC, get_elec);

    % Verbose diagnostics: print per-participant pair counts and expected
    verbose = true;

    region_fields = {'peripheral','central','central_peripheral','peripheral_peripheral','central_central'};
    region_labels = {'Peripheral','Central','Central-Peripheral','Peripheral-Peripheral','Central-Central'};

    window_defs = struct( ...
        'name', {'After_-15_0', 'After_0_-15'}, ...
        'row_min', {-15, 0}, ...
        'col_min', {0, -15});

    out_dir = fullfile(basefold, 'Results', 'Group_Analysis', 'Post_Point_Significance');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    summary_rows = [];
    row_idx = 0;

    for con = 1:numel(all_con)
        cond_name = char(all_con(con));
        fprintf('Analyzing post-point significance for condition %s\n', cond_name);

        n_windows = 1;
        if ~combine_windows
            n_windows = numel(window_defs);
        end
        region_summary = initialize_region_summary(numel(region_fields), n_windows);

        for p = 1:numel(participants)
            participant_name = char(participants(p));
            fprintf('  Loading %s (%d/%d)\n', participant_name, p, numel(participants));

            CoIData = max_get_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, cond_name, cutoff, timing);

            for r = 1:numel(region_fields)
                region_name = region_fields{r};
                if ~isfield(CoIData, region_name) || isempty(CoIData.(region_name).mask)
                    continue;
                end

                % Diagnostic: compare actual pairs found vs expected from EoI
                if verbose
                    if isfield(CoIData.(region_name), 'n_pairs')
                        actual_pairs = CoIData.(region_name).n_pairs;
                    else
                        actual_pairs = size(CoIData.(region_name).mask, 3);
                    end
                    expected_pairs = NaN;
                    % expected for single-electrode within-region fields
                    switch region_name
                        case 'peripheral'
                            if isfield(CoIData, 'peripheral_elecs')
                                expected_pairs = numel(CoIData.peripheral_elecs);
                            end
                        case 'central'
                            if isfield(CoIData, 'central_elecs')
                                expected_pairs = numel(CoIData.central_elecs);
                            end
                        case 'central_peripheral'
                            if isfield(CoIData, 'central_elecs') && isfield(CoIData, 'peripheral_elecs')
                                expected_pairs = numel(CoIData.central_elecs) * numel(CoIData.peripheral_elecs);
                            end
                        case 'peripheral_peripheral'
                            if isfield(CoIData, 'peripheral_elecs')
                                n = numel(CoIData.peripheral_elecs);
                                expected_pairs = n * (n - 1) / 2;
                            end
                        case 'central_central'
                            if isfield(CoIData, 'central_elecs')
                                n = numel(CoIData.central_elecs);
                                expected_pairs = n * (n - 1) / 2;
                            end
                    end
                    if ~isnan(expected_pairs) && actual_pairs ~= expected_pairs
                        fprintf('   [DEBUG] %s: region %s has %d pairs (expected %d)\n', participant_name, region_name, actual_pairs, expected_pairs);
                    end
                end

                [sig_counts, total_counts] = count_postpoint_significance(CoIData.(region_name).mask, timing, window_defs, frac_sig_threshold, combine_windows);
                region_summary(r).sig_counts = region_summary(r).sig_counts + sig_counts(:)';
                region_summary(r).total_counts = region_summary(r).total_counts + total_counts(:)';
            end
        end
        % Prepare window names (combined or separate)
        if combine_windows
            wn = {window_defs.name};
            combined_name = strjoin(wn, '_AND_');
            n_write_windows = 1;
        else
            combined_name = '';
            n_write_windows = numel(window_defs);
        end

        for r = 1:numel(region_fields)
            for w = 1:n_write_windows
                row_idx = row_idx + 1;
                sig_pairs = region_summary(r).sig_counts(w);
                total_pairs = region_summary(r).total_counts(w);
                pct_pairs = 100 * sig_pairs / max(total_pairs, 1);

                summary_rows(row_idx).Condition = string(cond_name); %#ok<AGROW>
                summary_rows(row_idx).Region = string(region_fields{r}); %#ok<AGROW>
                summary_rows(row_idx).Label = string(region_labels{r}); %#ok<AGROW>
                if combine_windows
                    summary_rows(row_idx).Window = string(combined_name); %#ok<AGROW>
                    summary_rows(row_idx).RowMinMs = NaN; %#ok<AGROW>
                    summary_rows(row_idx).ColMinMs = NaN; %#ok<AGROW>
                else
                    summary_rows(row_idx).Window = string(window_defs(w).name); %#ok<AGROW>
                    summary_rows(row_idx).RowMinMs = window_defs(w).row_min; %#ok<AGROW>
                    summary_rows(row_idx).ColMinMs = window_defs(w).col_min; %#ok<AGROW>
                end
                summary_rows(row_idx).SigPairs = sig_pairs; %#ok<AGROW>
                summary_rows(row_idx).TotalPairs = total_pairs; %#ok<AGROW>
                summary_rows(row_idx).PctPairs = pct_pairs; %#ok<AGROW>
            end
        end

        stats_this_cond = struct2table(summary_rows(strcmp({summary_rows.Condition}, cond_name)));
        csv_path = fullfile(out_dir, sprintf('postpoint_significance_%s.csv', cond_name));
        writetable(stats_this_cond, csv_path);
        fprintf('Saved summary table to %s\n', csv_path);

        plot_postpoint_summary(stats_this_cond, region_labels, window_defs, cond_name, out_dir, combine_windows, combined_name);
    end

    % Save combined results across all conditions
    if ~isempty(summary_rows)
        all_stats = struct2table(summary_rows);
        combined_csv = fullfile(out_dir, 'postpoint_significance_all_conditions.csv');
        writetable(all_stats, combined_csv);
        fprintf('Saved combined summary table to %s\n', combined_csv);

        % Also save as .mat for quick loading
        combined_mat = fullfile(out_dir, 'postpoint_significance_all_conditions.mat');
        save(combined_mat, 'all_stats');
        fprintf('Saved combined MAT to %s\n', combined_mat);
    end

    stats = struct2table(summary_rows);
end

function region_summary = initialize_region_summary(n_regions, n_windows)
    template = struct('sig_counts', zeros(1, n_windows), 'total_counts', zeros(1, n_windows));
    region_summary = repmat(template, n_regions, 1);
end

function [sig_counts, total_counts] = count_postpoint_significance(mask_stack, timing, window_defs, frac_sig_threshold, combine_windows)
    if ndims(mask_stack) == 2
        mask_stack = reshape(mask_stack, size(mask_stack,1), size(mask_stack,2), 1);
    end
    if combine_windows
        % Build union of all window masks
        combined_mask = false(length(timing), length(timing));
        for w = 1:numel(window_defs)
            wm = build_window_mask(timing, window_defs(w).row_min, window_defs(w).col_min);
            combined_mask = combined_mask | wm;
        end

        n_cells = nnz(combined_mask);
        sig_counts = 0;
        total_counts = 0;
        if n_cells > 0
            n_pairs = size(mask_stack, 3);
            total_counts = n_pairs;
            for p = 1:n_pairs
                pair_mask = mask_stack(:,:,p) > 0;
                frac_sig = nnz(pair_mask & combined_mask) / n_cells;
                if frac_sig > frac_sig_threshold
                    sig_counts = sig_counts + 1;
                end
            end
        end
        sig_counts = sig_counts(:)';
        total_counts = total_counts(:)';
    else
        n_windows = numel(window_defs);
        sig_counts = zeros(1, n_windows);
        total_counts = zeros(1, n_windows);

        for w = 1:n_windows
            window_mask = build_window_mask(timing, window_defs(w).row_min, window_defs(w).col_min);
            n_cells = nnz(window_mask);
            if n_cells == 0
                continue;
            end

            n_pairs = size(mask_stack, 3);
            total_counts(w) = n_pairs;

            for p = 1:n_pairs
                pair_mask = mask_stack(:,:,p) > 0;
                frac_sig = nnz(pair_mask & window_mask) / n_cells;
                if frac_sig > frac_sig_threshold
                    sig_counts(w) = sig_counts(w) + 1;
                end
            end
        end
    end
end

function window_mask = build_window_mask(timing, row_min, col_min)
    row_keep = timing >= row_min;
    col_keep = timing >= col_min;
    window_mask = row_keep(:) & col_keep(:)';
end

function plot_postpoint_summary(stats_this_cond, region_labels, window_defs, cond_name, out_dir, combine_windows, combined_name)
    if combine_windows
        n_tiles = 1;
    else
        n_tiles = numel(window_defs);
    end

    fig = figure('Color', 'w', 'Position', [100 100 1200 500]);
    tiledlayout(fig, 1, n_tiles, 'Padding', 'compact', 'TileSpacing', 'compact');

    for w = 1:n_tiles
        nexttile;
        if combine_windows
            rows = strcmp(stats_this_cond.Window, string(combined_name));
            title_str = combined_name;
        else
            rows = strcmp(stats_this_cond.Window, string(window_defs(w).name));
            title_str = window_defs(w).name;
        end
        data = stats_this_cond(rows, :);
        pct = data.PctPairs;
        bars = bar(pct, 'FaceColor', [0.35 0.55 0.80], 'EdgeColor', 'none'); %#ok<NASGU>
        ylim([0 100]);
        grid on;
        box off;
        xticks(1:numel(region_labels));
        xticklabels(region_labels);
        xtickangle(20);
        ylabel('Significant pair instances (%)');
        title(sprintf('%s', title_str), 'Interpreter', 'none');

        for r = 1:height(data)
            text(r, min(100, pct(r) + 4), sprintf('%d/%d', data.SigPairs(r), data.TotalPairs(r)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9);
        end
    end

    sgtitle(sprintf('Post-point CoI significance across flies (%s)', cond_name), 'Interpreter', 'none');
    saveas(fig, fullfile(out_dir, sprintf('postpoint_significance_%s.fig', cond_name)));
    exportgraphics(fig, fullfile(out_dir, sprintf('postpoint_significance_%s.png', cond_name)), 'Resolution', 300);
    try
        exportgraphics(fig, fullfile(out_dir, sprintf('postpoint_significance_%s.svg', cond_name)), 'ContentType', 'vector');
    catch
        print(fig, fullfile(out_dir, sprintf('postpoint_significance_%s.svg', cond_name)), '-dsvg');
    end
    close(fig);
end