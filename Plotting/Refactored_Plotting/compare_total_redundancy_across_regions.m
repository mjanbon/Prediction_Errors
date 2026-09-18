function compare_total_redundancy_across_conditions(participants, basefold, datatype, cutoff, times)
    % Compare total redundancy for sleep and wake across all regions

    % Define conditions and region labels
    activity_tags = {'wake', 'sleep'};
    region_fields = {'peripheral', 'central', 'central_peripheral', 'peripheral_peripheral', 'central_central'};
    region_labels = {'Peripheral', 'Central', 'Central-Peripheral', 'Peripheral-Peripheral', 'Central-Central'};
    n_regions = numel(region_fields);
    n_conditions = numel(activity_tags);

    % Preallocate results
    total_redundancy = zeros(n_regions, n_conditions);

    % Loop over conditions
    for c = 1:n_conditions
        activity_tag = activity_tags{c};
        % Use the first available condition in all_con (e.g., 'BSLEEP' or 'BWAKE')
        [~, ~, all_con, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = Max_get_param(0, 0);
        condition = char(all_con(1)); % Adjust if you want to loop over all_con

        % Get group-averaged region data
        CoI_group = get_group_region_data(participants, basefold, datatype, activity_tag, condition, cutoff, times);

        % For each region, sum the positive values (redundancy)
        for r = 1:n_regions
            region_name = region_fields{r};
            region_data = CoI_group.(region_name);
            if isfield(region_data, 'redundant') && ~isempty(region_data.redundant)
                total_redundancy(r, c) = sum(region_data.redundant(:), 'omitnan');
            else
                total_redundancy(r, c) = NaN;
            end
        end
    end

    % Plot results
    figure;
    bar(total_redundancy);
    set(gca, 'XTickLabel', region_labels, 'XTick', 1:n_regions);
    legend(activity_tags, 'Location', 'Best');
    ylabel('Total Redundancy (sum of positive CoI)');
    title('Total Redundancy for Sleep and Wake Across Regions');

    % Save figure
    save_dir = fullfile(basefold, 'Results', 'Group_Analysis');
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    saveas(gcf, fullfile(save_dir, 'Total_Redundancy_Sleep_Wake.png'));
    saveas(gcf, fullfile(save_dir, 'Total_Redundancy_Sleep_Wake.fig'));
end

function CoI_group = get_group_region_data(participants, basefold, datatype, activity_tag, condition, cutoff, times)
    % Aggregate region data across all flies for a given condition
    CoI_group = struct();
    region_fields = {'peripheral', 'central', 'central_peripheral', 'peripheral_peripheral', 'central_central'};
    for r = 1:numel(region_fields)
        CoI_group.(region_fields{r}) = initialize_combined_region_struct(times);
    end

    for subject = 1:numel(participants)
        participant_name = char(participants(subject));
        CoIData = max_get_plotting_CoI_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, times);

        % Accumulate for each region
        CoI_group.peripheral = accumulate_region_data(CoI_group.peripheral, CoIData.peripheral);
        CoI_group.central = accumulate_region_data(CoI_group.central, CoIData.central);
        CoI_group.central_peripheral = accumulate_region_data(CoI_group.central_peripheral, CoIData.central_peripheral);
        CoI_group.peripheral_peripheral = accumulate_region_data(CoI_group.peripheral_peripheral, CoIData.peripheral_peripheral);
        CoI_group.central_central = accumulate_region_data(CoI_group.central_central, CoIData.central_central);
    end

    % Average across flies
    for r = 1:numel(region_fields)
        CoI_group.(region_fields{r}) = average_region_data(CoI_group.(region_fields{r}), numel(participants));
    end
end