%% File: max_plot_main_all_flies_rf.m
% Group-level Co-I plotting driver that averages each region across flies

function max_plot_main_all_flies_rf(comparison_folder)
    set_default_plotting();

    if nargin < 1 || isempty(comparison_folder)
        comparison_folder = 'Comparisons';
    end

    %% Parameters
    USING_HPC = 0;
    get_elec = 0;
    cutoff = 10;
    timing = -25:1:74; % 100 timepoints
    times = 100;

    % Plot layout configuration
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
        'colormap_gamma', 0.65, ...
        'climits_mask', [0 1] ...
    );

    %% Load participant info
    [basefold, datatype, all_con, ~, ~, participants, ~, ~, activity_tag, ~, ~, ~, ~, ~, ~, ~, ~] = Max_get_param(USING_HPC, get_elec);

    plotted_fields = {'peripheral','central','central_peripheral','peripheral_peripheral','central_central'};
    plotted_labels = {'Peripheral','Central','Central-Peripheral','Peripheral-Peripheral','Central-Central'};

    %% Iterate over conditions and build group averages across flies
    for condition = 1:numel(all_con)
        cond_name = char(all_con(condition));
        fprintf('Averaging Co-I across %d flies for condition %s\n', numel(participants), cond_name);

        group_data = struct();
        for i = 1:numel(plotted_fields)
            group_data.(plotted_fields{i}) = initialize_group_region_struct(times);
        end

        for subject = 1:numel(participants)
            participant_name = char(participants(subject));
            fprintf('  Loading %s (%d/%d)\n', participant_name, subject, numel(participants));

            CoIData = max_get_plotting_CoI_rf(basefold, datatype, participant_name, activity_tag, cond_name, cutoff, times, comparison_folder);

            for i = 1:numel(plotted_fields)
                region_name = plotted_fields{i};
                if ~isfield(CoIData, region_name) || isempty(CoIData.(region_name).mean)
                    continue;
                end
                group_data.(region_name) = accumulate_group_region_data(group_data.(region_name), CoIData.(region_name));
            end
        end

        % Finalize all region data and compute a single global climits across regions
        valid_regions = false(1, numel(plotted_fields));
        max_vals = nan(1, numel(plotted_fields));
        for i = 1:numel(plotted_fields)
            region_name = plotted_fields{i};
            if group_data.(region_name).n_flies == 0
                fprintf('  Skipping %s: no valid flies found\n', plotted_labels{i});
                continue;
            end
            group_data.(region_name) = finalize_group_region_data(group_data.(region_name));
            vals = group_data.(region_name).mean(:);
            vals = vals(isfinite(vals));
            if ~isempty(vals)
                max_vals(i) = max(abs(vals));
                valid_regions(i) = true;
            end
        end

        if any(valid_regions)
            global_max = max(max_vals(valid_regions));
            global_climits = [-global_max, global_max];
        else
            global_climits = plot_cfg.climits;
        end

        % Plot each region using the shared global climits
        for i = 1:numel(plotted_fields)
            region_name = plotted_fields{i};
            if group_data.(region_name).n_flies == 0
                continue;
            end
            plot_cfg_curr = plot_cfg;
            plot_cfg_curr.climits = global_climits;
            fprintf('  Plotting %s with shared climits [%.6f %.6f]\n', plotted_labels{i}, global_climits(1), global_climits(2));
            plot_CoI_region(group_data.(region_name), basefold, 'All_Flies', activity_tag, cond_name, timing, plot_cfg_curr, plotted_labels{i});
        end
    end
end

function S = initialize_group_region_struct(times)
    n = numel(times);
    S = struct(...
        'mean_sum', zeros(n), ...
        'mask_sum', zeros(n), ...
        'mask_r_sum', zeros(n), ...
        'mask_s_sum', zeros(n), ...
        'mi1_all', [], ...
        'mi2_all', [], ...
        'n_pairs', 0, ...
        'n_flies', 0, ...
        'mean', [], ...
        'redundant', [], ...
        'synergetic', [], ...
        'mask', [], ...
        'mask_r', [], ...
        'mask_s', [], ...
        'mi1', [], ...
        'mi2', []);
end

function S_out = accumulate_group_region_data(S_out, S_in)
    if isempty(S_in.mean) || ~isfield(S_in, 'n') || S_in.n == 0
        return;
    end

    S_out.mean_sum = S_out.mean_sum + S_in.mean;
    S_out.mask_sum = S_out.mask_sum + double(S_in.mask);
    S_out.mask_r_sum = S_out.mask_r_sum + double(S_in.mask_r);
    S_out.mask_s_sum = S_out.mask_s_sum + double(S_in.mask_s);

    % Keep every electrode-pair MI trace so the plot error band reflects
    % variability across all pairs pooled across all flies.
    if ~isempty(S_in.mi1)
        S_out.mi1_all = [S_out.mi1_all, S_in.mi1]; %#ok<AGROW>
    end
    if ~isempty(S_in.mi2)
        S_out.mi2_all = [S_out.mi2_all, S_in.mi2]; %#ok<AGROW>
    end

    S_out.n_pairs = S_out.n_pairs + S_in.n;
    S_out.n_flies = S_out.n_flies + 1;
end

function S = finalize_group_region_data(S)
    if S.n_flies == 0
        return;
    end

    S.mean = S.mean_sum / S.n_flies;
    S.redundant = S.mean .* (S.mean > 0);
    S.synergetic = S.mean .* (S.mean < 0);

    S.mask = S.mask_sum;
    S.mask_r = S.mask_r_sum;
    S.mask_s = S.mask_s_sum;

    if ~isempty(S.mi1_all)
        S.mi1 = S.mi1_all;
    else
        S.mi1 = [];
    end
    if ~isempty(S.mi2_all)
        S.mi2 = S.mi2_all;
    else
        S.mi2 = [];
    end

    % Use the total number of contributing pair instances so the mask
    % panels show proportions across all flies and all available pairs.
    if S.n_pairs == 0
        S.n_pairs = 1;
    end
    S.n = S.n_pairs;
end

% compute_region_climits removed — global climits computed directly in caller