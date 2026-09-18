
%% File: max_plot_main.m
% Main script to load, iterate over subjects, and plot CoI analysis results

function max_plot_main()
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
        'climits_mask', [0 1] ...
    );

    %% Load participant info
    [basefold, datatype, all_con, ~, ~, participants, ~, ~, activity_tag, ~, ~, ~, ~, ~, ~, ~, ~] = Max_get_param(USING_HPC, get_elec);

    %% Iterate over all participants and conditions
    for subject = 1:numel(participants)
        participant_name = char(participants(subject));
        for condition = 1:numel(all_con)
            cond_name = char(all_con(condition));

            fprintf('Processing %s, condition %s\n', participant_name, cond_name);

            CoIData = max_get_plotting_CoI(basefold, datatype, participant_name, activity_tag, cond_name, cutoff, times);

            % Plot each region (within and cross-electrode types)
            if ~isempty(CoIData.peripheral.mean)
                plot_CoI_region(CoIData.peripheral, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg, 'Peripheral');
            end
            if ~isempty(CoIData.central.mean)
                plot_CoI_region(CoIData.central, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg, 'Central');
            end
            if ~isempty(CoIData.central_peripheral.mean)
                plot_CoI_region(CoIData.central_peripheral, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg, 'Central-Peripheral');
            end
            if ~isempty(CoIData.peripheral_peripheral.mean)
                plot_CoI_region(CoIData.peripheral_peripheral, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg, 'Peripheral-Peripheral');
            end
            if ~isempty(CoIData.central_central.mean)
                plot_CoI_region(CoIData.central_central, basefold, participant_name, activity_tag, cond_name, timing, plot_cfg, 'Central-Central');
            end
        end
    end
end
