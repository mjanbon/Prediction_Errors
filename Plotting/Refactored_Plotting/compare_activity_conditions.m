function results = compare_activity_conditions(fly_list, basefold, datatype, condition, cutoff, times, activity_cond1, activity_cond2, sig_threshold, dir_threshold, coi_source)
% COMPARE_ACTIVITY_CONDITIONS
% Uses CoI comparison outputs (diffCoI + sigMask) produced by
% max_main_comparison/max_Get_COI_comparison to summarize and plot
% activity-tag differences across flies.
%
% INPUTS:
%   fly_list        - Cell array of fly names to include
%   basefold        - Base folder for data
%   datatype        - Type of data being analyzed
%   condition       - Condition name (e.g., 'BSLEEP')
%   cutoff          - Electrode cutoff for central/peripheral division
%   times           - Time vector for the analysis
%   activity_cond1  - First activity condition (e.g., 'wake')
%   activity_cond2  - Second activity condition (e.g., 'sleep')
%   sig_threshold   - Threshold for high-consensus significance (default: 0.3)
%   dir_threshold   - Threshold for high directional consensus (default: 0.3)
%   coi_source      - 'pretrim' or 'original' comparison output (default: 'pretrim')
%
% OUTPUT:
%   results         - Struct containing mean diffCoI and consensus sigMask

    if nargin < 11 || isempty(coi_source)
        coi_source = 'pretrim';
    end
    if nargin < 10 || isempty(dir_threshold)
        dir_threshold = 0.3;
    end
    if nargin < 9 || isempty(sig_threshold)
        sig_threshold = 0.3;
    end
    if nargin < 8 || isempty(activity_cond2)
        activity_cond2 = 'sleep';
    end
    if nargin < 7 || isempty(activity_cond1)
        activity_cond1 = 'wake';
    end

    comparison_tag = [activity_cond1 '_' activity_cond2];
    regions = {'central', 'peripheral', 'central_peripheral', 'central_central', 'peripheral_peripheral'};

    fprintf('Comparing %s vs %s for condition %s using %s CoI\n', activity_cond1, activity_cond2, condition, coi_source);

    % Load EoI metadata. Pretrim comparison outputs have their own
    % pretrim-specific EoI file, so use that when coi_source asks for
    % pretrim data.
    EoI = load_eoi(basefold, datatype, coi_source);

    % Initialize aggregation containers
    results = struct();
    results.missing_eoi = {};
    results.missing_comparison_data = {};
    diff_stack = struct();
    sig_stack = struct();
    mi1_diff_stack = struct();
    mi2_diff_stack = struct();
    for r = 1:numel(regions)
        region = regions{r};
        diff_stack.(region) = {};
        sig_stack.(region) = {};
        mi1_diff_stack.(region) = {};
        mi2_diff_stack.(region) = {};
        results.(region).n_valid = 0;
    end

    % Iterate flies and aggregate region-level diffs
    for i = 1:numel(fly_list)
        pname = fly_list{i};

        elecs = get_eoi_for_fly(EoI, pname, comparison_tag, condition);
        if isempty(elecs)
            fprintf('No EoI entry for %s (%s), skipping\n', pname, comparison_tag);
            results.missing_eoi(end+1, :) = {pname, comparison_tag, condition};
            continue;
        end

        CoIComp = load_comparison_data(basefold, pname, activity_cond1, activity_cond2, condition, coi_source);
        if isempty(CoIComp)
            warning('No %s comparison data loaded for %s %s. Skipping this fly.', coi_source, pname, condition);
            results.missing_comparison_data(end+1, :) = {pname, condition, coi_source};
            continue;
        end

        [central_elecs, peripheral_elecs] = split_elecs(elecs, cutoff);

        region_defs = struct();
        region_defs.central = {central_elecs, central_elecs, true};
        region_defs.peripheral = {peripheral_elecs, peripheral_elecs, true};
        region_defs.central_peripheral = {central_elecs, peripheral_elecs, false};
        region_defs.central_central = {central_elecs, central_elecs, true};
        region_defs.peripheral_peripheral = {peripheral_elecs, peripheral_elecs, true};

        for r = 1:numel(regions)
            region = regions{r};
            def = region_defs.(region);
            elecsA = def{1};
            elecsB = def{2};
            skip_duplicates = def{3};

            if isempty(elecsA) || isempty(elecsB)
                continue;
            end

            [diff_mean, sig_mean, mi1_diff, mi2_diff, n_pairs] = process_region_pairs(CoIComp, pname, condition, elecsA, elecsB, skip_duplicates, times);
            if n_pairs == 0
                continue;
            end

            diff_stack.(region){end+1} = diff_mean;
            sig_stack.(region){end+1} = sig_mean;
            mi1_diff_stack.(region){end+1} = mi1_diff;
            mi2_diff_stack.(region){end+1} = mi2_diff;
            results.(region).n_valid = results.(region).n_valid + 1;
        end
    end

    % Final aggregation across flies
    for r = 1:numel(regions)
        region = regions{r};
        if isempty(diff_stack.(region))
            results.(region).diff_mean = [];
            results.(region).sig_mask = [];
            results.(region).mi1_diff = [];
            results.(region).mi2_diff = [];
            continue;
        end

        diff_mat = cat(3, diff_stack.(region){:});
        sig_mat = cat(3, sig_stack.(region){:});
        mi1_diff_mat = cat(3, mi1_diff_stack.(region){:});
        mi2_diff_mat = cat(3, mi2_diff_stack.(region){:});

        results.(region).diff_mean = mean(diff_mat, 3, 'omitnan');
        
        % Calculate proportion of electrode pairs (across all flies) that show significance
        % sig_mat is [time x time x n_flies*n_pairs_per_fly]
        % We simply average across all instances to get proportion
        mean_sig = mean(sig_mat, 3, 'omitnan');
        results.(region).sig_mask = mean_sig;  % Proportion ranges from -1 to +1
        
        results.(region).mi1_diff = mean(mi1_diff_mat, 3, 'omitnan');
        results.(region).mi2_diff = mean(mi2_diff_mat, 3, 'omitnan');
    end

    if ~isempty(results.missing_eoi)
        fprintf('Missing EoI entries for %d fly/condition rows:\n', size(results.missing_eoi, 1));
        disp(results.missing_eoi);
    end
    if ~isempty(results.missing_comparison_data)
        fprintf('Missing %s comparison data for %d fly/condition rows:\n', coi_source, size(results.missing_comparison_data, 1));
        disp(results.missing_comparison_data);
    end

    % Plot results
    plot_comparison_diff(results, times, activity_cond1, activity_cond2, condition, basefold, sig_threshold, dir_threshold);
end

function EoI = load_eoi(basefold, datatype, coi_source)
    if nargin < 3
        coi_source = 'original';
    end

    if any(strcmpi(coi_source, {'pretrim', 'coi_pretrim', 'comparison_pretrim'}))
        eoi_file_candidates = { ...
            ['EoI_data_' datatype '_pretrim.mat'], ...
            ['EoI_data_' datatype '.mat']};
    else
        eoi_file_candidates = { ...
            ['EoI_data_' datatype '.mat'], ...
            ['EoI_data_' datatype '_pretrim.mat']};
    end

    eoi_dir_candidates = { ...
        fullfile(basefold, 'DataEoI'), ...
        fullfile(basefold, 'DataEOI')};

    eoi_path = '';
    for dir_ind = 1:numel(eoi_dir_candidates)
        for file_ind = 1:numel(eoi_file_candidates)
            candidate_path = fullfile(eoi_dir_candidates{dir_ind}, eoi_file_candidates{file_ind});
            if exist(candidate_path, 'file')
                eoi_path = candidate_path;
                break;
            end
        end
        if ~isempty(eoi_path)
            break;
        end
    end

    if isempty(eoi_path)
        error('EoI file not found for datatype %s and coi_source %s', datatype, coi_source);
    end

    fprintf('Loading EoI from %s\n', eoi_path);
    tmp = load(eoi_path, 'EoI');
    EoI = tmp.EoI;
end

function elecs = get_eoi_for_fly(EoI, pname, comparison_tag, condition)
    elecs = {};
    if isfield(EoI, pname) && isfield(EoI.(pname), comparison_tag) && isfield(EoI.(pname).(comparison_tag), condition)
        elecs = EoI.(pname).(comparison_tag).(condition);
        return;
    end

    % Fallback: try inverted tag order
    parts = split(comparison_tag, '_');
    if numel(parts) == 2
        alt_tag = [parts{2} '_' parts{1}];
        if isfield(EoI, pname) && isfield(EoI.(pname), alt_tag) && isfield(EoI.(pname).(alt_tag), condition)
            elecs = EoI.(pname).(alt_tag).(condition);
        end
    end
end

function [central_elecs, peripheral_elecs] = split_elecs(all_elecs, cutoff)
    if isempty(all_elecs)
        central_elecs = {};
        peripheral_elecs = {};
        return;
    end
    if isscalar(cutoff)
        central_elecs = all_elecs(1:cutoff-1);
        peripheral_elecs = all_elecs(cutoff:end);
    else
        error('cutoff must be a scalar for this comparison function');
    end
end

function CoIComp = load_comparison_data(basefold, pname, activity_cond1, activity_cond2, condition, coi_source)
    CoIComp = struct();
    loaded_any = false;

    switch lower(coi_source)
        case {'pretrim', 'coi_pretrim', 'comparison_pretrim'}
            comparison_dir = 'Comparisons_pretrim';
            variable_names = {'CoI_pretrim', 'CoI', 'MI_stat'};
        case {'original', 'coi', 'raw'}
            comparison_dir = 'Comparisons';
            variable_names = {'CoI', 'MI_stat', 'CoI_pretrim'};
        otherwise
            error('Unknown coi_source "%s". Use "pretrim" or "original".', coi_source);
    end

    comp_folder = [pname '_' activity_cond1 '_' activity_cond2 '_' condition];
    comp_dir = find_comparison_dir(basefold, comparison_dir, comp_folder);
    if ~exist(comp_dir, 'dir')
        CoIComp = [];
        return;
    end

    files = dir(fullfile(comp_dir, '*.mat'));
    if isempty(files)
        CoIComp = [];
        return;
    end

    for i = 1:numel(files)
        file_path = fullfile(comp_dir, files(i).name);
        file_base = erase(files(i).name, '.mat');

        data = load(file_path);
        src = extract_comparison_source(data, file_base, variable_names);
        if isempty(src) || ~isfield(src, pname) || ~isfield(src.(pname), condition)
            continue;
        end

        % Build electrode-pair identifier
        identifier = erase(file_base, [pname '_']);
        identifier = erase(identifier, 'permuted_');

        cond_data = src.(pname).(condition);
        if isfield(cond_data, 'diffCoI') && isfield(cond_data.diffCoI, identifier)
            CoIComp.(pname).(condition).diffCoI.(identifier) = cond_data.diffCoI.(identifier);
            loaded_any = true;
        end
        if isfield(cond_data, 'sigMask') && isfield(cond_data.sigMask, identifier)
            CoIComp.(pname).(condition).sigMask.(identifier) = cond_data.sigMask.(identifier);
        end
        if isfield(cond_data, 'MI1_diff') && isfield(cond_data.MI1_diff, identifier)
            CoIComp.(pname).(condition).mi1_diff.(identifier) = cond_data.MI1_diff.(identifier);
        elseif isfield(cond_data, 'mi1_diff') && isfield(cond_data.mi1_diff, identifier)
            CoIComp.(pname).(condition).mi1_diff.(identifier) = cond_data.mi1_diff.(identifier);
        end
        if isfield(cond_data, 'MI2_diff') && isfield(cond_data.MI2_diff, identifier)
            CoIComp.(pname).(condition).mi2_diff.(identifier) = cond_data.MI2_diff.(identifier);
        elseif isfield(cond_data, 'mi2_diff') && isfield(cond_data.mi2_diff, identifier)
            CoIComp.(pname).(condition).mi2_diff.(identifier) = cond_data.mi2_diff.(identifier);
        end
    end

    if ~loaded_any
        warning('Found %s comparison files for %s %s, but none contained matching diffCoI fields.', coi_source, pname, condition);
        CoIComp = [];
    end
end

function comp_dir = find_comparison_dir(basefold, comparison_dir, comp_folder)
    candidates = {
        fullfile(basefold, 'Drosophila_CoI', comparison_dir, comp_folder)
        fullfile(basefold, 'Data', 'Drosophila_CoI', comparison_dir, comp_folder)
        fullfile(fileparts(basefold), 'Drosophila_CoI', comparison_dir, comp_folder)
    };

    comp_dir = candidates{1};
    for c = 1:numel(candidates)
        if exist(candidates{c}, 'dir')
            comp_dir = candidates{c};
            return;
        end
    end
end

function src = extract_comparison_source(data, file_base, variable_names)
    src = [];
    for v = 1:numel(variable_names)
        var_name = variable_names{v};
        if ~isfield(data, var_name)
            continue;
        end

        candidate = data.(var_name);
        if strcmp(var_name, 'MI_stat') && isfield(candidate, 'CoI')
            src = candidate.CoI;
            return;
        end
        if isfield(candidate, file_base) && isfield(candidate.(file_base), 'CoI')
            src = candidate.(file_base).CoI;
            return;
        end
        if isfield(candidate, 'CoI')
            src = candidate.CoI;
            return;
        end
    end
end

function [diff_mean, sig_mean, mi1_diff, mi2_diff, n_pairs] = process_region_pairs(CoIComp, pname, cond, elecsA, elecsB, skip_duplicates, times)
    n = length(times);
    FFi_all = [];
    sig_all = [];
    mi1_all = [];
    mi2_all = [];
    elec_pairs = {};

    if isempty(CoIComp) || ~isfield(CoIComp, pname) || ~isfield(CoIComp.(pname), cond) || ...
            ~isfield(CoIComp.(pname).(cond), 'diffCoI')
        diff_mean = [];
        sig_mean = [];
        mi1_diff = [];
        mi2_diff = [];
        n_pairs = 0;
        return;
    end

    cond_data = CoIComp.(pname).(cond);

    for i = 1:numel(elecsA)
        for j = 1:numel(elecsB)
            E1 = char(elecsA(i));
            E2 = char(elecsB(j));
            if strcmp(E1, E2) && ~isequal(elecsA, elecsB)
                continue;
            end
            pair_name = [E1 '_' E2];
            pair_inv = [E2 '_' E1];

            if skip_duplicates && ismember(pair_inv, elec_pairs)
                continue;
            end
            elec_pairs{end+1} = pair_name;

            if myIsField(cond_data.diffCoI, pair_name)
                FFi_all = cat(3, FFi_all, cond_data.diffCoI.(pair_name));
                if myIsField(cond_data, 'sigMask') && myIsField(cond_data.sigMask, pair_name)
                    sig_all = cat(3, sig_all, cond_data.sigMask.(pair_name));
                else
                    sig_all = cat(3, sig_all, zeros(n, n));
                end
                if myIsField(cond_data, 'mi1_diff') && myIsField(cond_data.mi1_diff, pair_name)
                    mi1_all = cat(2, mi1_all, cond_data.mi1_diff.(pair_name)(:));
                end
                if myIsField(cond_data, 'mi2_diff') && myIsField(cond_data.mi2_diff, pair_name)
                    mi2_all = cat(2, mi2_all, cond_data.mi2_diff.(pair_name)(:));
                end
            end
        end
    end

    n_pairs = size(FFi_all, 3);
    if n_pairs == 0
        diff_mean = [];
        sig_mean = [];
        mi1_diff = [];
        mi2_diff = [];
        return;
    end

    diff_mean = mean(FFi_all, 3, 'omitnan');
    sig_mean = mean(sig_all, 3, 'omitnan');
    if isempty(mi1_all)
        mi1_diff = nan(n, 1);
    else
        mi1_diff = mean(mi1_all, 2, 'omitnan');
    end
    if isempty(mi2_all)
        mi2_diff = nan(n, 1);
    else
        mi2_diff = mean(mi2_all, 2, 'omitnan');
    end
end

function plot_comparison_diff(results, times, cond1, cond2, data_condition, basefold, sig_threshold, dir_threshold)
    % Calculate global color limits for CoI difference
    all_diffs = [];
    region_labels = {'central', 'peripheral', 'central_peripheral', 'central_central', 'peripheral_peripheral'};
    for r = 1:numel(region_labels)
        region = region_labels{r};
        if isfield(results, region) && ~isempty(results.(region).diff_mean)
            all_diffs = [all_diffs; results.(region).diff_mean(:)];
        end
    end
    if isempty(all_diffs)
        warning('No comparison data available to plot for %s vs %s - %s.', cond1, cond2, data_condition);
        return;
    end
    diff_clim = max(abs([min(all_diffs), max(all_diffs)]));
    if diff_clim == 0
        diff_clim = 1;
    end

    % Create separate figure for each region (2x3 layout)
    region_display = {'Central', 'Peripheral', 'Central-Peripheral', 'Central-Central', 'Peripheral-Peripheral'};

    save_dir = fullfile(basefold, 'Results', 'Condition_Comparison');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    for r = 1:numel(region_labels)
        region = region_labels{r};
        region_disp = region_display{r};

        if ~isfield(results, region) || isempty(results.(region).diff_mean)
            fprintf('Skipping %s (no data)\n', region);
            continue;
        end

        fig = figure('Position', [100 100 1800 600]);
        tiledlayout(2, 3, 'TileSpacing', 'Compact');
        sgtitle(sprintf('%s vs %s - %s (%s)', cond1, cond2, region_disp, data_condition));

        sig_mask = results.(region).sig_mask;
        diff_data = results.(region).diff_mean;

        % Row 1, Position 1: Main difference CoI matrix
        ax1 = nexttile(1);
        contourf(ax1, times, times, diff_data, 50, 'linecolor', 'none');
        colormap(ax1, redblue(256));
        shading(ax1, 'flat');
        caxis(ax1, [-diff_clim, diff_clim]);
        colorbar(ax1);
        daspect(ax1, [1 1 1]);
        xlabel(ax1, 'Time (ms)');
        ylabel(ax1, 'Time (ms)');
        title(ax1, sprintf('CoI Difference (%s - %s)', cond2, cond1));

        % Row 1, Position 2: Significance mask
        ax_sig = nexttile(2);
        imagesc(ax_sig, times, times, sig_mask);
        set(ax_sig, 'YDir', 'normal');
        colormap(ax_sig, 'gray');
        caxis(ax_sig, [-1, 1]);
        colorbar(ax_sig);
        daspect(ax_sig, [1 1 1]);
        xlabel(ax_sig, 'Time (ms)');
        ylabel(ax_sig, 'Time (ms)');
        title(ax_sig, 'Proportion Significant');
        
        % Add text showing overall proportion significant
            % total_sig = sum(abs(sig_mask(:)) > 0) / numel(sig_mask);
            % text(ax_sig, min(times), max(times)*0.95, sprintf('%.1f%% sig', total_sig*100), ...
            %     'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'k');

        % Row 1, Position 3: Thresholded significance (|proportion| > sig_threshold)
        ax_thresh = nexttile(3);
        sig_thresh = sig_mask;
        sig_thresh(abs(sig_mask) <= sig_threshold) = 0;  % Zero out points below threshold
        imagesc(ax_thresh, times, times, sig_thresh);
        set(ax_thresh, 'YDir', 'normal');
        colormap(ax_thresh, flipud(gray));
        caxis(ax_thresh, [-1, 1]);
        colorbar(ax_thresh);
        daspect(ax_thresh, [1 1 1]);
        xlabel(ax_thresh, 'Time (ms)');
        ylabel(ax_thresh, 'Time (ms)');
        title(ax_thresh, sprintf('High Consensus (|prop| > %.2f)', sig_threshold));
        
        % Add text showing proportion of high-consensus points
            % high_consensus = sum(abs(sig_mask(:)) > sig_threshold) / numel(sig_mask);
            % text(ax_thresh, min(times), max(times)*0.95, sprintf('%.1f%% high', high_consensus*100), ...
            %     'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'k');

        % Row 2, Position 1: Directional consensus (proportion agreeing with mean CoI direction)
        ax_dir = nexttile(4);
        % For each time point, check if mean CoI diff is positive or negative
        % Then count proportion of sig_mask values that match that direction
        dir_consensus = zeros(size(diff_data));
        for ti = 1:size(diff_data, 1)
            for tj = 1:size(diff_data, 2)
                if diff_data(ti, tj) > 0
                    % Mean CoI shows cond2 > cond1, count proportion where sig_mask = +1
                    % sig_mask ranges from -1 to +1, where +1 means all pairs show cond2>cond1
                    % We want proportion that are positive (agreeing with direction)
                    dir_consensus(ti, tj) = max(0, sig_mask(ti, tj));
                elseif diff_data(ti, tj) < 0
                    % Mean CoI shows cond1 > cond2, count proportion where sig_mask = -1
                    % We want proportion that are negative (agreeing with direction)
                    dir_consensus(ti, tj) = max(0, -sig_mask(ti, tj));
                else
                    % No difference in mean, set to 0
                    dir_consensus(ti, tj) = 0;
                end
            end
        end
        imagesc(ax_dir, times, times, dir_consensus);
        set(ax_dir, 'YDir', 'normal');
        colormap(ax_dir, flipud(gray));
        caxis(ax_dir, [0, 1]);
        colorbar(ax_dir);
        daspect(ax_dir, [1 1 1]);
        xlabel(ax_dir, 'Time (ms)');
        ylabel(ax_dir, 'Time (ms)');
        title(ax_dir, 'Directional Consensus');
        
        % Add text showing overall directional consensus
            % mean_dir_consensus = mean(dir_consensus(:), 'omitnan');
            % text(ax_dir, min(times), max(times)*0.95, sprintf('%.1f%% agree', mean_dir_consensus*100), ...
            %     'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'k');

        % Row 2, Position 2: High Directional Consensus (dir_consensus > dir_threshold)
        ax_high_dir = nexttile(5);
        high_dir_consensus = dir_consensus;
        high_dir_consensus(dir_consensus <= dir_threshold) = 0;  % Zero out points below threshold
        imagesc(ax_high_dir, times, times, high_dir_consensus);
        set(ax_high_dir, 'YDir', 'normal');
        colormap(ax_high_dir, flipud(gray));
        caxis(ax_high_dir, [0, 1]);
        colorbar(ax_high_dir);
        daspect(ax_high_dir, [1 1 1]);
        xlabel(ax_high_dir, 'Time (ms)');
        ylabel(ax_high_dir, 'Time (ms)');
        title(ax_high_dir, sprintf('High Directional Consensus (> %.2f)', dir_threshold));
        
        % Add text showing proportion of high directional consensus points
            % high_dir_prop = sum(dir_consensus(:) > dir_threshold) / numel(dir_consensus);
            % text(ax_high_dir, min(times), max(times)*0.95, sprintf('%.1f%% high', high_dir_prop*100), ...
            %     'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'k');

        % Save figure
        fig_file = fullfile(save_dir, sprintf('%s_vs_%s_%s_diff.fig', cond1, cond2, region));
        png_file = fullfile(save_dir, sprintf('%s_vs_%s_%s_diff.png', cond1, cond2, region));
        saveas(fig, fig_file, 'fig');
        saveas(fig, png_file, 'png');
        fprintf('Saved %s\n', region);
        close(fig);
    end
end

function tf = myIsField(s, f)
    tf = isfield(s, f);
end
