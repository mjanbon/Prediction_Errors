%% File: max_plot_main_rf.m
% Main script to load, iterate over subjects, and plot CoI analysis results

function max_plot_main_rf()
    set_default_plotting();

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
    

    %% Iterate over all participants and conditions
    plotted_fields = {'peripheral','central','central_peripheral','peripheral_peripheral','central_central'};
    for subject = 1:numel(participants)
        participant_name = char(participants(subject));

        % Preload all conditions and compute one subject-level Co-I max so
        % every region plot for this subject shares the same color limits.
        subject_max_abs = NaN;
        CoIData_cache = cell(numel(all_con), 1);
        cond_name_cache = strings(numel(all_con), 1);
        for condition = 1:numel(all_con)
            cond_name = char(all_con(condition));
            cond_name_cache(condition) = string(cond_name);
            CoIData = max_get_plotting_CoI_rf(basefold, datatype, participant_name, activity_tag, cond_name, cutoff, times);
            CoIData_cache{condition} = CoIData;
            local_max_abs = compute_global_coi_max_abs(CoIData, plotted_fields);
            if ~isnan(local_max_abs)
                if isnan(subject_max_abs)
                    subject_max_abs = local_max_abs;
                else
                    subject_max_abs = max(subject_max_abs, local_max_abs);
                end
            end
        end

        if isnan(subject_max_abs) || subject_max_abs <= 0
            subject_climits = plot_cfg.climits;
        else
            subject_climits = [-subject_max_abs, subject_max_abs];
        end
        fprintf('Using subject-level Co-I clim for %s: [%.6f %.6f]\n', participant_name, subject_climits(1), subject_climits(2));

        for condition = 1:numel(all_con)
            cond_name = char(cond_name_cache(condition));
            CoIData = CoIData_cache{condition};

            fprintf('Processing %s, condition %s\n', participant_name, cond_name);

            % Use one shared Co-I color scale across all regions for this
            % subject, so color saturation is directly comparable.
            plot_cfg_curr = plot_cfg;
            plot_cfg_curr.climits = subject_climits;
            fprintf('Using shared Co-I clim for %s/%s: [%.6f %.6f]\n', participant_name, cond_name, plot_cfg_curr.climits(1), plot_cfg_curr.climits(2));

            % Plot each region (within and cross-electrode types)
            if ~isempty(CoIData.peripheral.mean)
                plot_CoI_region(CoIData.peripheral, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg_curr, 'Peripheral');
            end
            if ~isempty(CoIData.central.mean)
                plot_CoI_region(CoIData.central, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg_curr, 'Central');
            end
            if ~isempty(CoIData.central_peripheral.mean)
                plot_CoI_region(CoIData.central_peripheral, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg_curr, 'Central-Peripheral');
            end
            if ~isempty(CoIData.peripheral_peripheral.mean)
                plot_CoI_region(CoIData.peripheral_peripheral, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg_curr, 'Peripheral-Peripheral');
            end
            if ~isempty(CoIData.central_central.mean)
                plot_CoI_region(CoIData.central_central, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg_curr, 'Central-Central');
            end
            
            % Add full matrix analysis by distance
            try
                fprintf('Running full matrix analysis by distance...\n');
                analyze_full_matrix_synergy(basefold, datatype, participant_name, activity_tag, cond_name, timing);
            catch e
                fprintf('Error in full matrix analysis: %s\n', e.message);
                disp(getReport(e, 'extended'));
            end
        end
    end
end

function max_abs = compute_global_coi_max_abs(CoIData, region_fields)
% Compute max absolute Co-I across available region mean matrices.

max_abs = NaN;
for i = 1:numel(region_fields)
    f = region_fields{i};
    if ~isfield(CoIData, f)
        continue;
    end
    region = CoIData.(f);
    if ~isfield(region, 'mean') || isempty(region.mean)
        continue;
    end
    vals = region.mean(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    local_max = max(abs(vals));
    if isnan(max_abs)
        max_abs = local_max;
    else
        max_abs = max(max_abs, local_max);
    end
end
end
