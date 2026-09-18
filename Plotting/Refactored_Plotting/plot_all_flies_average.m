%% Parameters
    USING_HPC = 0;
    get_elec = 0;
    cutoff = 10;
    times = -25:1:74; % 100 timepoints
    activity_tags = {'wake', 'sleep'};
    compare_conditions = 'False'
    single_condition = 'False'
    co_i_by_distance = 'False'
    co_i_by_electrode = 'True'
    plot_separate_regions = 'False' % New flag for the separate region plots
    total_coi_across_conditions = 'False'
    compare_total_coi_across_4_conditions = 'False'
    do_compare_region_patterns = 'False'
    coi_source = 'pretrim'; % Options: 'pretrim' or 'original'
    
    % times = 99;

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
        'climits', [-0.03 0.03], ...
        'climits_mask', [0 1] ...
    );

    %% Load participant info
    [basefold, datatype, all_con, ~, ~, participants, ~, ~, activity_tag, ~, ~, ~, ~, ~, ~, ~, ~] = Max_get_param(USING_HPC, get_elec);
    condition = char(all_con(1))
    if strcmp(single_condition,'True')==1
        plot_group_average_CoI(participants, basefold, datatype, activity_tag, condition, cutoff, times)
    elseif strcmp(do_compare_region_patterns,'True')==1
        % R and SSIM between different processing regions
        compare_region_patterns(participants, basefold, datatype, activity_tag, condition, cutoff, times);
    elseif strcmp(compare_conditions,'True')==1
        % Example 1: Compare sleep vs. wake (default)
        compare_activity_conditions(participants, basefold, datatype, condition, cutoff, times, 'wake', 'beginning_sleep', [], [], coi_source);
    elseif strcmp(co_i_by_distance,'True')==1
        analyze_synergy_redundancy_across_flies(participants, basefold, datatype, activity_tag, condition, times);
    elseif strcmp(co_i_by_electrode,'True')==1
        analyze_synergy_redundancy_by_electrode(participants, basefold, datatype, activity_tag, condition);
    elseif strcmp(plot_separate_regions,'True')==1
        % New function to plot separate regions with MI
        plot_all_regions_with_MI(participants, basefold, datatype, activity_tag, condition, cutoff, times, plot_cfg);
    elseif strcmp(total_coi_across_conditions,'True')==1
        stats_all = compare_total_coi_across_conditions(participants, basefold, datatype, cutoff, times, activity_tags);
        stats = stats_all.redundancy;
    elseif strcmp(compare_total_coi_across_4_conditions,'True')==1
        stats_all = compare_total_coi_four_activity_tags(participants, basefold, datatype, cutoff, times)
    end

function plot_all_regions_with_MI(participants, basefold, datatype, activity_tag, condition, cutoff, times, plot_cfg)
    % Plot separate figures for each region type with both CoI and MI, averaged across all flies
    % Following the format used in plot_CoI_region.m
    
    % Initialize data containers for each region type
    peripheral_data = initialize_combined_region_struct(times);
    central_data = initialize_combined_region_struct(times);
    central_peripheral_data = initialize_combined_region_struct(times);
    peripheral_peripheral_data = initialize_combined_region_struct(times);
    central_central_data = initialize_combined_region_struct(times);
    
    % Initialize counters for total pairs per region
    total_pairs = struct(...
        'peripheral', 0, ...
        'central', 0, ...
        'central_peripheral', 0, ...
        'peripheral_peripheral', 0, ...
        'central_central', 0);
    
    % Load and combine data from all flies
    fprintf('Loading data from all flies...\n');
    for subject = 1:numel(participants)
        participant_name = char(participants(subject));
        fprintf('Processing %s, condition %s\n', participant_name, condition);
        
        % Get CoI data for this fly
        CoIData = max_get_plotting_CoI_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, times);
        
        % Store MI data for each fly to calculate SED later
        if ~isempty(CoIData.peripheral.mi1) && CoIData.peripheral.n > 0
            peripheral_data.mi1_all(:, subject) = mean(CoIData.peripheral.mi1, 2);
            peripheral_data.mi2_all(:, subject) = mean(CoIData.peripheral.mi2, 2);
            total_pairs.peripheral = total_pairs.peripheral + CoIData.peripheral.n;
        end
        
        if ~isempty(CoIData.central.mi1) && CoIData.central.n > 0
            central_data.mi1_all(:, subject) = mean(CoIData.central.mi1, 2);
            central_data.mi2_all(:, subject) = mean(CoIData.central.mi2, 2);
            total_pairs.central = total_pairs.central + CoIData.central.n;
        end
        
        if ~isempty(CoIData.central_peripheral.mi1) && CoIData.central_peripheral.n > 0
            central_peripheral_data.mi1_all(:, subject) = mean(CoIData.central_peripheral.mi1, 2);
            central_peripheral_data.mi2_all(:, subject) = mean(CoIData.central_peripheral.mi2, 2);
            total_pairs.central_peripheral = total_pairs.central_peripheral + CoIData.central_peripheral.n;
        end
        
        if ~isempty(CoIData.peripheral_peripheral.mi1) && CoIData.peripheral_peripheral.n > 0
            peripheral_peripheral_data.mi1_all(:, subject) = mean(CoIData.peripheral_peripheral.mi1, 2);
            peripheral_peripheral_data.mi2_all(:, subject) = mean(CoIData.peripheral_peripheral.mi2, 2);
            total_pairs.peripheral_peripheral = total_pairs.peripheral_peripheral + CoIData.peripheral_peripheral.n;
        end
        
        if ~isempty(CoIData.central_central.mi1) && CoIData.central_central.n > 0
            central_central_data.mi1_all(:, subject) = mean(CoIData.central_central.mi1, 2);
            central_central_data.mi2_all(:, subject) = mean(CoIData.central_central.mi2, 2);
            total_pairs.central_central = total_pairs.central_central + CoIData.central_central.n;
        end
        
        % Accumulate data for each region type
        peripheral_data = accumulate_region_data(peripheral_data, CoIData.peripheral);
        central_data = accumulate_region_data(central_data, CoIData.central);
        central_peripheral_data = accumulate_region_data(central_peripheral_data, CoIData.central_peripheral);
        peripheral_peripheral_data = accumulate_region_data(peripheral_peripheral_data, CoIData.peripheral_peripheral);
        central_central_data = accumulate_region_data(central_central_data, CoIData.central_central);
    end
    
    % Calculate averages across all flies
    peripheral_data = average_region_data(peripheral_data, numel(participants));
    central_data = average_region_data(central_data, numel(participants));
    central_peripheral_data = average_region_data(central_peripheral_data, numel(participants));
    peripheral_peripheral_data = average_region_data(peripheral_peripheral_data, numel(participants));
    central_central_data = average_region_data(central_central_data, numel(participants));
    
    % Set total electrode pairs
    peripheral_data.n = total_pairs.peripheral;
    central_data.n = total_pairs.central;
    central_peripheral_data.n = total_pairs.central_peripheral;
    peripheral_peripheral_data.n = total_pairs.peripheral_peripheral;
    central_central_data.n = total_pairs.central_central;
    
    % Get common timing vector for x-axis
    timing = -25:1:74;
    
    % Create plot output directory if it doesn't exist
    plot_dir = fullfile(basefold, 'Plots', 'Group_Average_Regions');
    if ~exist(plot_dir, 'dir')
        mkdir(plot_dir);
    end
    
    % Plot each region type using the style from plot_CoI_region.m
    fprintf('Creating plots for each region type...\n');
    
    plot_regional_average(peripheral_data, plot_dir, 'All_Flies', activity_tag, condition, timing, plot_cfg, 'Peripheral');
    plot_regional_average(central_data, plot_dir, 'All_Flies', activity_tag, condition, timing, plot_cfg, 'Central');
    plot_regional_average(central_peripheral_data, plot_dir, 'All_Flies', activity_tag, condition, timing, plot_cfg, 'Central-Peripheral');
    plot_regional_average(peripheral_peripheral_data, plot_dir, 'All_Flies', activity_tag, condition, timing, plot_cfg, 'Peripheral-Peripheral');
    plot_regional_average(central_central_data, plot_dir, 'All_Flies', activity_tag, condition, timing, plot_cfg, 'Central-Central');
    
    fprintf('All plots created successfully in %s\n', plot_dir);
end

function S = initialize_combined_region_struct(times)
    n = length(times);
    S = struct(...
        'FFi_all', [], ...
        'mask_sum', zeros(n, n), ...
        'mask_count', zeros(n, n), ...
        'mask_r_sum', zeros(n, n), ...
        'mask_r_count', zeros(n, n), ...
        'mask_s_sum', zeros(n, n), ...
        'mask_s_count', zeros(n, n), ...
        'mi1_sum', zeros(1, n), ...
        'mi2_sum', zeros(1, n), ...
        'mi1_all', nan(n, 8), ...  % Store MI for each fly (max 8)
        'mi2_all', nan(n, 8), ...  % Store MI for each fly (max 8)
        'mean_sum', zeros(n), ...
        'redundant_sum', zeros(n), ...
        'synergetic_sum', zeros(n), ...
        'n_elecs', 0, ...
        'n_flies', 0);
end

function S_out = accumulate_region_data(S_out, S_in)
    % Accumulate data from one fly to the combined structure
    if isempty(S_in.mean) || isempty(S_in.mi1) || S_in.n == 0
        return;  % Skip empty data
    end

    % Sum the matrices
    S_out.mean_sum = S_out.mean_sum + S_in.mean;
    S_out.redundant_sum = S_out.redundant_sum + S_in.redundant;
    S_out.synergetic_sum = S_out.synergetic_sum + S_in.synergetic;

    % --- Correct mask accumulation for all pairs and all flies ---
    if ~isempty(S_in.mask)
        if ndims(S_in.mask) == 3
            % S_in.mask is [T x T x n_pairs]
            for pairIdx = 1:size(S_in.mask, 3)
                mask_this = double(S_in.mask(:,:,pairIdx));
                valid = ~isnan(mask_this);
                S_out.mask_sum(valid) = S_out.mask_sum(valid) + mask_this(valid);
                S_out.mask_count(valid) = S_out.mask_count(valid) + 1;
            end
        elseif ismatrix(S_in.mask)
            % S_in.mask is [T x T] and ALREADY averaged over pairs for this fly
            % Just sum over flies, not pairs!
            S_out.mask_sum = S_out.mask_sum + double(S_in.mask);
            S_out.mask_count = S_out.mask_count + 1;
        end
    end

    if ~isempty(S_in.mask_r)
        if ndims(S_in.mask_r) == 3
            S_out.mask_r_sum = S_out.mask_r_sum + sum(double(S_in.mask_r), 3);
            S_out.mask_r_count = S_out.mask_r_count + size(S_in.mask_r, 3);
        elseif ismatrix(S_in.mask_r)
            S_out.mask_r_sum = S_out.mask_r_sum + double(S_in.mask_r);
            S_out.mask_r_count = S_out.mask_r_count + 1;
        end
    end

    if ~isempty(S_in.mask_s)
        if ndims(S_in.mask_s) == 3
            S_out.mask_s_sum = S_out.mask_s_sum + sum(double(S_in.mask_s), 3);
            S_out.mask_s_count = S_out.mask_s_count + size(S_in.mask_s, 3);
        elseif ismatrix(S_in.mask_s)
            S_out.mask_s_sum = S_out.mask_s_sum + double(S_in.mask_s);
            S_out.mask_s_count = S_out.mask_s_count + 1;
        end
    end

    % Sum MI data (using mean for each fly)
    S_out.mi1_sum = S_out.mi1_sum + mean(S_in.mi1, 2)';
    S_out.mi2_sum = S_out.mi2_sum + mean(S_in.mi2, 2)';

    % Count electrodes and flies
    S_out.n_elecs = S_out.n_elecs + S_in.n;
    S_out.n_flies = S_out.n_flies + 1;
end

function S = average_region_data(S, n_flies)
    % Calculate average values across flies
    if S.n_flies == 0
        return;  % Avoid division by zero
    end

    % Average the matrices
    S.mean = S.mean_sum / S.n_flies;
    S.redundant = S.redundant_sum / S.n_flies;
    S.synergetic = S.synergetic_sum / S.n_flies;

    % Calculate mask averages - proportion of all electrode pairs that were significant
    mask_valid = S.mask_count > 0;
    S.mask = zeros(size(S.mask_sum));
    S.mask(mask_valid) = S.mask_sum(mask_valid) ./ S.mask_count(mask_valid);

    mask_r_valid = S.mask_r_count > 0;
    S.mask_r = zeros(size(S.mask_r_sum));
    S.mask_r(mask_r_valid) = S.mask_r_sum(mask_r_valid) ./ S.mask_r_count(mask_r_valid);

    mask_s_valid = S.mask_s_count > 0;
    S.mask_s = zeros(size(S.mask_s_sum));
    S.mask_s(mask_s_valid) = S.mask_s_sum(mask_s_valid) ./ S.mask_s_count(mask_s_valid);

    % Average MI data
    S.mi1 = S.mi1_sum / S.n_flies;
    S.mi2 = S.mi2_sum / S.n_flies;
end

function plot_regional_average(region, plot_dir, pname, tag, cond, timing, cfg, label)
    % Create a tiled plot of CoI and MI results for a given region, averaged across flies
    % This follows the exact format used in plot_CoI_region.m

    % Skip if no data available
    if isempty(region.mean) || ~any(region.mean(:))
        fprintf('No data available for %s region.\n', label);
        return;
    end

    % Create figure with tiled layout - match exactly plot_CoI_region.m
    fig = figure('Position', [100 100 800 600], 'Renderer', 'painters');
    tiledlayout(4, 3, 'TileSpacing', 'Compact');
    sgtitle(sprintf('Group Average %s %s %s CoI', pname, tag, label));  % Global figure title

    plot_coi_sign_key(cfg);
    plot_label_tile(7, 'Redundancy only', [0.80 0.10 0.10]);
    plot_label_tile(10, 'Synergy only', [0.05 0.20 0.85]);

    % Main CoI matrix - position 5 (middle)
    ax1 = nexttile(5);
    contourf(ax1, timing, timing, region.mean, 50, 'linecolor', 'none');
    clim(cfg.climits); colormap(ax1, redblue(256)); shading(ax1, 'flat');
    cb = colorbar; daspect(ax1, [1 1 1]);
    xlim(ax1, cfg.xlimits); ylim(ax1, cfg.ylimits_CoI);
    xlabel('Time (ms)'); ylabel('Time (ms)');
    title(ax1, [label ' CoI']);
    set(ax1, 'XTick', cfg.xticks_CoI, 'XTickLabel', cfg.x_labels, ...
        'YTick', cfg.yticks_CoI, 'YTickLabel', cfg.y_labels_CoI);
    cb.Label.String = 'Co-I (bits)';
    cb.Label.FontSize = 11;
    cb.Label.Rotation = 270;
    cb.Label.VerticalAlignment = 'bottom';

    % Redundancy & Synergy CoI
    ax_red = nexttile(8);
    contourf(ax_red, timing, timing, region.redundant, 50, 'linecolor', 'none');
    clim(cfg.climits); colormap(ax_red, redblue(256)); shading(ax_red, 'flat');
    cb_red = colorbar; daspect(ax_red, [1 1 1]);
    cb_red.Label.String = 'Co-I (bits)';
    xlim(ax_red, cfg.xlimits); ylim(ax_red, cfg.ylimits_CoI);
    xlabel('Time'); ylabel('Time'); 
    title(ax_red, [label ' Redundant']);
    set(ax_red, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);

    ax_syn = nexttile(11);
    contourf(ax_syn, timing, timing, region.synergetic, 50, 'linecolor', 'none');
    clim(cfg.climits); colormap(ax_syn, redblue(256)); shading(ax_syn, 'flat');
    cb_syn = colorbar; daspect(ax_syn, [1 1 1]);
    cb_syn.Label.String = 'Co-I (bits)';
    xlim(ax_syn, cfg.xlimits); ylim(ax_syn, cfg.ylimits_CoI);
    xlabel('Time'); ylabel('Time'); 
    title(ax_syn, [label ' Synergetic']);
    set(ax_syn, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);

    % Masks with n annotation - match exactly what's in plot_CoI_region.m
    ax_mask = nexttile(6);
    imagesc(ax_mask, timing, timing, region.mask); 
    set(ax_mask, 'YDir', 'normal'); daspect(ax_mask, [1 1 1]);
    clim([0 1]); colormap(ax_mask, flipud(bone(256)));
    cb_mask = colorbar(ax_mask);
    cb_mask.Label.String = 'Proportion significant';
    xlim(ax_mask, cfg.xlimits); ylim(ax_mask, cfg.ylimits_CoI);
    xlabel('Time'); ylabel('Time');
    title(ax_mask, sprintf('Significance mask (n = %d)', region.n));
    set(ax_mask, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);

    ax_rmask = nexttile(9);
    imagesc(ax_rmask, timing, timing, region.mask_r); 
    set(ax_rmask, 'YDir', 'normal'); daspect(ax_rmask, [1 1 1]);
    clim([0 1]); colormap(ax_rmask, flipud(bone(256)));
    cb_rmask = colorbar(ax_rmask);
    cb_rmask.Label.String = 'Proportion significant';
    xlim(ax_rmask, cfg.xlimits); ylim(ax_rmask, cfg.ylimits_CoI);
    xlabel('Time'); ylabel('Time');
    title(ax_rmask, 'Redundancy mask');
    set(ax_rmask, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);

    ax_smask = nexttile(12);
    imagesc(ax_smask, timing, timing, region.mask_s); 
    set(ax_smask, 'YDir', 'normal'); daspect(ax_smask, [1 1 1]);
    clim([0 1]); colormap(ax_smask, flipud(bone(256)));
    cb_smask = colorbar(ax_smask);
    cb_smask.Label.String = 'Proportion significant';
    xlim(ax_smask, cfg.xlimits); ylim(ax_smask, cfg.ylimits_CoI);
    xlabel('Time'); ylabel('Time');
    title(ax_smask, 'Synergy mask');
    set(ax_smask, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);

    % MI1 and MI2 curves - Use stdshade-like visualization
    ax2 = nexttile(2);
    plot(ax2, timing(:), region.mi1(:), 'g-', 'LineWidth', 0.5);
    hold(ax2, 'on');
    valid_flies = sum(~isnan(region.mi1_all(1,:)));
    if valid_flies > 0
        mi1_std = nanstd(region.mi1_all, 0, 2);
        mi1_sed = mi1_std / sqrt(valid_flies);
        fill([timing(:); flipud(timing(:))], ...
             [region.mi1(:) + mi1_sed; flipud(region.mi1(:) - mi1_sed)], ...
             'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    hold(ax2, 'off');
    title(ax2, 'MI1');
    xlim(ax2, cfg.xlimits); ylim(ax2, cfg.ylimit_MI); 
    xlabel(ax2, 'Time'); ylabel(ax2, 'MI (bits)');
    set(ax2, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_MI);
    grid(ax2, 'on');

    % MI2 plot - robust to shape and always displays
    ax4 = nexttile(4);
    mi2 = region.mi2(:); % force column
    timing = timing(:);  % force column
    if length(mi2) ~= length(timing)
        warning('MI2 and timing length mismatch!');
        mi2 = interp1(linspace(1,length(mi2),length(mi2)), mi2, linspace(1,length(mi2),length(timing)), 'linear', 'extrap')';
    end
    plot(ax4, timing, mi2, 'r-', 'LineWidth', 0.5);
    hold(ax4, 'on');
    if size(region.mi2_all,2) > 0
        valid_cols = any(~isnan(region.mi2_all), 1);
        mi2_std = nanstd(region.mi2_all(:,valid_cols), 0, 2);
        mi2_sed = mi2_std / sqrt(sum(valid_cols));
        fill([timing; flipud(timing)], ...
             [mi2 + mi2_sed; flipud(mi2 - mi2_sed)], ...
             'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    hold(ax4, 'off');
    title(ax4, 'MI2');
    xlim(ax4, cfg.xlimits); ylim(ax4, cfg.ylimit_MI); 
    xlabel(ax4, 'Time'); ylabel(ax4, 'MI (bits)');
    set(ax4, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_MI);
    grid(ax4, 'on');

    % Debugging information for MI2
    disp(['MI2 size: ', mat2str(size(region.mi2))]);
    disp(['MI2 min/max: ', num2str(min(region.mi2(:))), ' / ', num2str(max(region.mi2(:)))]);
    disp(['Any NaN in MI2? ', num2str(any(isnan(region.mi2(:))))]);
    disp(['MI2_all size: ', mat2str(size(region.mi2_all))]);
    disp(['Any NaN in MI2_all? ', num2str(any(isnan(region.mi2_all(:))))]);

    % Save the figure
    filename = sprintf('Group_Average_%s_%s_%s.png', label, tag, cond);
    saveas(fig, fullfile(plot_dir, filename));
    saveas(fig, fullfile(plot_dir, strrep(filename, '.png', '.fig')));
    close(fig);

    fprintf('Created plot for %s region.\n', label);
end

function plot_CoI_tile(data, timing, cfg, tile_pos, title_text, span)
    if nargin < 6 || isempty(span)
        ax = nexttile(tile_pos);
    else
        ax = nexttile(tile_pos, span);
    end
    contourf(ax, timing, timing, data, 50, 'linecolor', 'none');
    clim(cfg.climits); colormap(ax, redblue(256)); shading(ax, 'flat');
    cb = colorbar; daspect(ax, [1 1 1]);
    cb.Label.String = 'Co-I (bits)';
    xlim(ax, cfg.xlimits); ylim(ax, cfg.ylimits_CoI);
    xlabel('Time'); ylabel('Time'); title(ax, title_text);
    set(ax, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);
end

function plot_coi_sign_key(cfg)
    ax = nexttile(1);
    axis(ax, 'off');
    clim(ax, cfg.climits);
    colormap(ax, redblue(256));
    cb = colorbar(ax, 'eastoutside');
    cb.Ticks = [cfg.climits(1), cfg.climits(2)];
    cb.TickLabels = {'Synergistic', 'Redundant'};
    cb.Label.String = 'Co-I (bits)';
end

function plot_label_tile(tile_pos, label_text, color_value)
    ax = nexttile(tile_pos);
    axis(ax, 'off');
    text(ax, 0.5, 0.5, label_text, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', color_value, ...
        'FontWeight', 'bold', ...
        'FontSize', 12);
end

function plot_mask_tile(data, timing, cfg, tile_pos, title_text, n, span)
    if nargin < 7 || isempty(span)
        ax = nexttile(tile_pos);
    else
        ax = nexttile(tile_pos, span);
    end
    imagesc(ax, timing, timing, data); set(ax, 'YDir', 'normal'); daspect(ax, [1 1 1]);
    clim([0 1]); colormap(ax, flipud(bone(256)));
    cb = colorbar;
    cb.Label.String = 'Proportion significant';
    xlim(ax, cfg.xlimits); ylim(ax, cfg.ylimits_CoI);
    xlabel('Time'); ylabel('Time');
    title(ax, sprintf('%s (n = %d)', title_text, n));
    set(ax, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);
end
