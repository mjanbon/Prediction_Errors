function plot_all_flies_average_2(participants, basefold, datatype, activity_tag, condition, cutoff, timing, plot_cfg)
% Group-level CoI/MI/mask plots using RAW per-pair data from all flies

% Determine region groups dynamically from the RAW CoIData structure.
% We call the RAW loader for the first participant that returns a non-empty
% structure and detect within-region fields (e.g. 'central','intermediate','peripheral').
regions = {};
region_labels = {};
found = false;
for subject = 1:numel(participants)
    participant_name = char(participants(subject));
    CoIData = max_get_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, timing);
    % look for candidate within-region fields (no underscore, non-meta)
    fn = fieldnames(CoIData);
    cand = {};
    for f = 1:numel(fn)
        fname = fn{f};
        % skip helper/meta fields
        if any(strcmp(fname, {'central_elecs','peripheral_elecs'}))
            continue;
        end
        % consider fields that do not contain an underscore and that have FFi_all
        if isempty(strfind(fname,'_')) && isfield(CoIData.(fname),'FFi_all')
            cand{end+1} = fname; %#ok<AGROW>
        end
    end
    if ~isempty(cand)
        % prefer ordered canonical list if present
        canonical = {'central','intermediate','peripheral'};
        regions = intersect(canonical, cand, 'stable');
        % if canonical intersection empty, use cand as-is
        if isempty(regions)
            regions = cand;
        end
        % build pretty labels
        region_labels = cellfun(@(s) regexprep(s,'_','-'), regions, 'UniformOutput', false);
        region_labels = cellfun(@(s) regexprep(s,'(^.|-.)', '${upper($0)}', 'once'), region_labels, 'UniformOutput', false);
        found = true;
        break;
    end
end
if ~found
    error('Could not detect any within-region fields from max_get_plotting_CoI_raw_rf outputs.');
end

% Initialize group data structures
n_regions = numel(regions);
for r = 1:n_regions
    group.(regions{r}).FFi_all = [];
    group.(regions{r}).mask_all = [];
    group.(regions{r}).mask_r_all = [];
    group.(regions{r}).mask_s_all = [];
    group.(regions{r}).mi1_all = [];
    group.(regions{r}).mi2_all = [];
    group.(regions{r}).n_pairs = 0;
end

% Aggregate all flies, all pairs using the detected region fields
contributors = {}; % participants that contributed any data
for subject = 1:numel(participants)
    participant_name = char(participants(subject));
    % Use the RAW version!
    CoIData = max_get_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, timing);
    contributed = false;

    for r = 1:n_regions
        reg = regions{r};
        if ~isfield(CoIData, reg)
            continue;
        end
        fly = CoIData.(reg);

        % Skip empty
        if isempty(fly) || ~isfield(fly,'FFi_all') || isempty(fly.FFi_all)
            continue;
        end

        contributed = true;

        % --- CoI matrices: concatenate all pairs from all flies ---
        group.(reg).FFi_all = cat(3, group.(reg).FFi_all, fly.FFi_all);

        % --- Masks: concatenate all pairs from all flies ---
        group.(reg).mask_all   = cat(3, group.(reg).mask_all, fly.mask);
        group.(reg).mask_r_all = cat(3, group.(reg).mask_r_all, fly.mask_r);
        group.(reg).mask_s_all = cat(3, group.(reg).mask_s_all, fly.mask_s);

        % --- MI1/MI2: concatenate all pairs from all flies ---
        if isfield(fly,'mi1'), group.(reg).mi1_all = [group.(reg).mi1_all, fly.mi1]; end
        if isfield(fly,'mi2'), group.(reg).mi2_all = [group.(reg).mi2_all, fly.mi2]; end

        group.(reg).n_pairs = group.(reg).n_pairs + size(fly.FFi_all,3);
    end

    if contributed
        contributors{end+1} = participant_name; %#ok<AGROW>
    end
end

% Compute group averages and plot
for r = 1:n_regions
    reg = regions{r};
    label = region_labels{r};
    if isempty(group.(reg).FFi_all), continue; end

    region.mean = mean(group.(reg).FFi_all, 3);
    region.redundant = region.mean .* (region.mean > 0);
    region.synergetic = region.mean .* (region.mean < 0);

    % --- Masks: total significant pair counts at each (i,j) across ALL pairs ---
    % We'll sum across the 3rd dimension (each slice is one electrode-pair)
    % and later divide by region.n (total number of pairs) when plotting
    % so the displayed value is the fraction of all pairs significant.
    region.mask   = sum(group.(reg).mask_all, 3, 'omitnan');
    region.mask_r = sum(group.(reg).mask_r_all, 3, 'omitnan');
    region.mask_s = sum(group.(reg).mask_s_all, 3, 'omitnan');
    region.n = group.(reg).n_pairs;

    % --- MI1/MI2: average across all pairs and flies ---
    region.mi1_all = group.(reg).mi1_all;
    region.mi2_all = group.(reg).mi2_all;
    region.mi1 = nanmean(group.(reg).mi1_all, 2);
    region.mi2 = nanmean(group.(reg).mi2_all, 2);

    % --- Plot using the same layout as plot_CoI_region ---
    plot_CoI_region_group(region, basefold, 'All_Flies', activity_tag, condition, timing, plot_cfg, label, n_regions, numel(contributors));
end

end

function plot_CoI_region_group(region, basefold, pname, tag, cond, timing, cfg, label, n_groups, n_flies)
% Tiled plot for group average, matching plot_CoI_region layout

% Set figure title and MI label/color logic based on region label
if contains(lower(label), 'central-peripheral')
    fig_title = 'Central-peripheral interactions Co-I';
    mi1_label = 'Peripheral MI';
    mi1_color = [0.8500 0.3250 0.0980]; % orange
    mi2_label = 'Central MI';
    mi2_color = [0.4660 0.6740 0.1880]; % green
elseif contains(lower(label), 'peripheral-peripheral')
    fig_title = 'Peripheral-peripheral interactions Co-I';
    mi1_label = 'Peripheral MI';
    mi1_color = [0.8500 0.3250 0.0980]; % orange
    mi2_label = 'Peripheral MI';
    mi2_color = [0.8500 0.3250 0.0980]; % orange
elseif contains(lower(label), 'central-central')
    fig_title = 'Central-central interactions Co-I';
    mi1_label = 'Central MI';
    mi1_color = [0.4660 0.6740 0.1880]; % green
    mi2_label = 'Central MI';
    mi2_color = [0.4660 0.6740 0.1880]; % green
elseif contains(lower(label), 'peripheral')
    fig_title = 'Peripheral regions Co-I';
    mi1_label = 'Peripheral MI';
    mi1_color = [0.8500 0.3250 0.0980]; % orange
    mi2_label = 'Peripheral MI';
    mi2_color = [0.8500 0.3250 0.0980]; % orange
elseif contains(lower(label), 'central')
    fig_title = 'Central regions Co-I';
    mi1_label = 'Central MI';
    mi1_color = [0.4660 0.6740 0.1880]; % green
    mi2_label = 'Central MI';
    mi2_color = [0.4660 0.6740 0.1880]; % green
else
    fig_title = [label ' Co-I'];
    mi1_label = 'MI1';
    mi1_color = [0 0 0];
    mi2_label = 'MI2';
    mi2_color = [0 0 0];
end

fig = figure('Position', [100 100 900 700], 'Renderer', 'painters');
tiledlayout(4,3,'TileSpacing','Compact');
sgtitle(fig_title);

plot_coi_sign_key(cfg);
plot_label_tile(7, 'Redundancy only', [0.80 0.10 0.10]);
plot_label_tile(10, 'Synergy only', [0.05 0.20 0.85]);

% MI1 (tile 2, top left)
ax2 = nexttile(2);
mi1 = region.mi1(:);
mi1_all = region.mi1_all;
if isempty(mi1_all), mi1_all = mi1; end
valid_cols = any(~isnan(mi1_all), 1);
mi1_std = nanstd(mi1_all(:,valid_cols), 0, 2);
mi1_sed = mi1_std / sqrt(sum(valid_cols));
plot(ax2, timing(:), mi1, '-', 'Color', mi1_color, 'LineWidth', 0.5); hold(ax2, 'on');
fill([timing(:); flipud(timing(:))], ...
     [mi1 + mi1_sed; flipud(mi1 - mi1_sed)], ...
     mi1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold(ax2, 'off');
title(ax2, mi1_label, 'Color', mi1_color, 'FontWeight', 'bold');
xlim(ax2, cfg.xlimits); ylim(ax2, cfg.ylimit_MI);
set(ax2, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_MI);
xlabel(ax2, 'Time (ms)'); ylabel(ax2, 'MI (bits)');
daspect(ax2, [1 1 1]);
axis(ax2, 'square');

% MI2 (tile 4, top right)
ax4 = nexttile(4);
mi2 = region.mi2(:);
if length(mi2) ~= length(timing)
    warning('MI2 and timing length mismatch! Interpolating MI2.');
    mi2 = interp1(linspace(1,length(mi2),length(mi2)), mi2, linspace(1,length(mi2),length(timing)), 'linear', 'extrap')';
end
plot(ax4, timing(:), mi2, '-', 'Color', mi2_color, 'LineWidth', 0.5); hold(ax4, 'on');
mi2_all = region.mi2_all;
if isempty(mi2_all), mi2_all = mi2; end
valid_cols = any(~isnan(mi2_all), 1);
if any(valid_cols)
    mi2_std = nanstd(mi2_all(:,valid_cols), 0, 2);
    mi2_sed = mi2_std / sqrt(sum(valid_cols));
    fill([timing(:); flipud(timing(:))], ...
         [mi2 + mi2_sed; flipud(mi2 - mi2_sed)], ...
         mi2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
hold(ax4, 'off');
view(ax4, -90, 90);
title(ax4, mi2_label, 'Color', mi2_color, 'FontWeight', 'bold');
xlim(ax4, cfg.xlimits); ylim(ax4, cfg.ylimit_MI);
set(ax4, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_MI);
xlabel(ax4, 'Time (ms)'); ylabel(ax4, 'MI (bits)');
daspect(ax4, [1 1 1]);
axis(ax4, 'square');

% Main CoI matrix (tile 5, center)
ax1 = nexttile(5);
contourf(ax1, timing, timing, region.mean, 50, 'linecolor', 'none');
clim(cfg.climits); colormap(ax1, redblue(256)); shading(ax1, 'flat');
cb = colorbar(ax1);
daspect(ax1, [1 1 1]);
xlim(ax1, cfg.xlimits); ylim(ax1, cfg.ylimits_CoI);
xlabel(ax1, 'Time (ms)'); ylabel(ax1, 'Time (ms)');
title(ax1, [label ' CoI']);
set(ax1,'XTick',cfg.xticks_CoI, 'XTickLabel', cfg.x_labels, ...
    'YTick', cfg.yticks_CoI, 'YTickLabel', cfg.y_labels_CoI);
cb.Label.String = 'Co-I (bits)';

% Redundancy & Synergy CoI (match layout from plot_CoI_region)
plot_CoI_tile(region.redundant, timing, cfg, 8, [label ' Redundant']);
plot_CoI_tile(region.synergetic, timing, cfg, 11, [label ' Synergetic']);

% Masks: all / redundancy / synergy (normalized by n)
plot_mask_tile(region.mask ./ max(1, region.n), timing, cfg, 6, 'All mask', region.n);
plot_mask_tile(region.mask_r ./ max(1, region.n), timing, cfg, 9, 'Redundancy mask', region.n);
plot_mask_tile(region.mask_s ./ max(1, region.n), timing, cfg, 12, 'Synergy mask', region.n);

% Save — include number of groups and contributing flies in filename
out_dir = fullfile(basefold, 'Plots', 'Group_Average_Regions');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if nargin < 9, n_groups = NaN; end
if nargin < 10, n_flies = NaN; end
filename = sprintf('Group_Average_g%d_f%d_%s_%s_%s', n_groups, n_flies, label, tag, cond);
saveas(fig, fullfile(out_dir, [filename '.fig']), 'fig');
saveas(fig, fullfile(out_dir, [filename '.jpg']), 'jpeg');
print(fig, fullfile(out_dir, [filename '.svg']), '-dsvg');
close(fig);
end

function plot_CoI_tile(data, timing, cfg, tile_pos, title_text)
    ax = nexttile(tile_pos);
    contourf(ax, timing, timing, data, 50, 'linecolor', 'none');
    clim(cfg.climits); colormap(ax, redblue(256)); shading(ax, 'flat');
    cb = colorbar; daspect(ax, [1 1 1]);
    cb.Label.String = 'Co-I (bits)';
    xlim(ax, cfg.xlimits); ylim(ax, cfg.ylimits_CoI);
    xlabel('Time (ms)'); ylabel(ax, 'Time (ms)'); title(ax, title_text);
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

function plot_mask_tile(data, timing, cfg, tile_pos, title_text, n)
    ax = nexttile(tile_pos);
    imagesc(ax, timing, timing, data); set(ax,'YDir','normal'); daspect(ax, [1 1 1]);
    clim(cfg.climits_mask); colormap(ax, flipud(bone(256)));
    cb = colorbar;
    cb.Label.String = 'Proportion significant';
    xlim(ax, cfg.xlimits); ylim(ax, cfg.ylimits_CoI);
    xlabel('Time (ms)'); ylabel(ax, 'Time (ms)');
    title(ax, sprintf('%s (n = %d)', title_text, n));
    set(ax, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);
end
