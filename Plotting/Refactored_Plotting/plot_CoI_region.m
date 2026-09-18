
function plot_CoI_region(region, basefold, pname, tag, cond, timing, cfg, label)
% PLOT_COI_REGION
% Creates a tiled plot of CoI and MI results for a given electrode group/region.

fig = figure('Position', [100 100 760 660], 'Renderer', 'painters');
tl = tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
sg = sgtitle(sprintf('%s %s %s CoI', pname, tag, label));  % Global figure title
set(sg, 'Interpreter', 'none');

if isfield(cfg, 'colormap_gamma') && ~isempty(cfg.colormap_gamma)
    cmap_gamma = cfg.colormap_gamma;
else
    cmap_gamma = 0.65;
end
left_column_shift = get_left_column_shift(cfg);

plot_coi_sign_key(tl, cfg, cmap_gamma);

% Main CoI matrix
ax1 = nexttile(tl, 5);
contourf(ax1, timing, timing, region.mean, 50, 'linecolor', 'none');
clim(cfg.climits); colormap(ax1, enhanced_redblue(256, cmap_gamma)); shading(ax1, 'flat');
xcb = colorbar(ax1);
set_standard_colorbar_ticks(xcb, cfg.climits, 5);
xcb.Label.String = 'Co-I (bits)';
tighten_colorbar(xcb);
daspect(ax1, [1 1 1]);
xlim(ax1, cfg.xlimits); ylim(ax1, cfg.ylimits_CoI);
xlabel('Time (ms)'); ylabel('Time (ms)');
title(ax1, [label ' CoI']);
lift_title(ax1);
set(ax1,'XTick',cfg.xticks_CoI, 'XTickLabel', cfg.x_labels, ...
    'YTick', cfg.yticks_CoI, 'YTickLabel', cfg.y_labels_CoI);

% Redundancy & Synergy CoI
plot_CoI_tile(tl, region.redundant, timing, cfg, 8, [label ' Redundant'], cmap_gamma);
plot_CoI_tile(tl, region.synergetic, timing, cfg, 11, [label ' Synergistic'], cmap_gamma);

% Masks with n annotation
plot_mask_tile(tl, region.mask ./ region.n, timing, cfg, 6, 'All mask', region.n);
plot_mask_tile(tl, region.mask_r ./ region.n, timing, cfg, 9, 'Redundancy mask', region.n);
plot_mask_tile(tl, region.mask_s ./ region.n, timing, cfg, 12, 'Synergy mask', region.n);

% MI1 and MI2 curves
ax2 = nexttile(tl, 2);
[mi1_col, mi2_col] = get_mi_colors(label);
mi1_title = get_mi_title(label, 1);
plot_mi_band(ax2, timing, region.mi1, mi1_col);
th = title(ax2, mi1_title); set(th,'Interpreter','none');
xlabel(ax2, 'Time (ms)'); ylabel(ax2, 'MI (bits)');
xlim(ax2, cfg.xlimits); ylim(ax2, cfg.ylimit_MI);
lift_title(ax2);
set(ax2, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_MI);
pbaspect(ax2, [1 1 1]);

ax4 = shifted_tile_axes(fig, tl, 4, left_column_shift);
mi2_title = get_mi_title(label, 2);
plot_mi_band(ax4, timing, region.mi2, mi2_col);
view(ax4, -90, 90); th2 = title(ax4, mi2_title); set(th2,'Interpreter','none');
xlabel(ax4, 'Time (ms)'); ylabel(ax4, 'MI (bits)');
xlim(ax4, cfg.xlimits); ylim(ax4, cfg.ylimit_MI);
lift_title(ax4);
set(ax4, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_MI);
pbaspect(ax4, [1 1 1]);

% Create the shifted text panels after the tiled axes/colorbars have settled.
drawnow;
plot_label_tile(fig, tl, 7, 'Redundancy only', [0.80 0.10 0.10], left_column_shift);
plot_label_tile(fig, tl, 10, 'Synergy only', [0.05 0.20 0.85], left_column_shift);

% Keep MATLAB's tiledlayout in control of positions so MI and CoI axes stay aligned.
drawnow;

% Save
save_plot(basefold, pname, tag, label, cond);
end

function plot_CoI_tile(tl, data, timing, cfg, tile_pos, title_text, cmap_gamma)
    ax = nexttile(tl, tile_pos);
    contourf(ax, timing, timing, data, 50, 'linecolor', 'none');
    clim(cfg.climits); colormap(ax, enhanced_redblue(256, cmap_gamma)); shading(ax, 'flat');
    xcb = colorbar(ax);
    set_standard_colorbar_ticks(xcb, cfg.climits, 5);
    xcb.Label.String = 'Co-I (bits)';
    tighten_colorbar(xcb);
    daspect(ax, [1 1 1]);
    xlim(ax, cfg.xlimits); ylim(ax, cfg.ylimits_CoI);
    xlabel('Time (ms)'); ylabel('Time (ms)'); title(ax, title_text);
    lift_title(ax);
    set(ax, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);
end


function plot_mask_tile(tl, data, timing, cfg, tile_pos, title_text, n)
    ax = nexttile(tl, tile_pos);
    imagesc(ax, timing, timing, data); set(ax,'YDir','normal'); daspect(ax, [1 1 1]);
    clim(cfg.climits_mask); colormap(ax, flipud(bone(256)));
    xcb = colorbar(ax);
    set_standard_colorbar_ticks(xcb, cfg.climits_mask, 3);
    xcb.Label.String = 'Proportion significant';
    tighten_colorbar(xcb);
    xlim(ax, cfg.xlimits); ylim(ax, cfg.ylimits_CoI);
    xlabel('Time (ms)'); ylabel('Time (ms)');
    title(ax, sprintf('%s (n = %d)', title_text, n));
    lift_title(ax);
    set(ax, 'XTick', cfg.xticks_CoI, 'YTick', cfg.yticks_CoI);
end

function plot_mi_band(ax, timing, values, line_color)
    if isempty(values)
        axes(ax); %#ok<LAXES>
        return;
    end

    if isvector(values)
        values = values(:);
    end

    if size(values, 1) ~= numel(timing) && size(values, 2) == numel(timing)
        values = values.';
    end

    mean_vec = mean(values, 2, 'omitnan');
    n_vec = sum(isfinite(values), 2);
    sem_vec = nanstd(values, 0, 2);
    sem_vec(n_vec > 0) = sem_vec(n_vec > 0) ./ sqrt(n_vec(n_vec > 0));
    sem_vec(n_vec == 0) = NaN;

    hold(ax, 'on');
    fill(ax, [timing(:); flipud(timing(:))], ...
        [mean_vec - sem_vec; flipud(mean_vec + sem_vec)], ...
        line_color, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
    plot(ax, timing(:), mean_vec, 'Color', line_color, 'LineWidth', 0.5);
    hold(ax, 'off');
end

function plot_coi_sign_key(tl, cfg, cmap_gamma)
    ax = nexttile(tl, 1);
    axis(ax, 'off');
    clim(ax, cfg.climits);
    colormap(ax, enhanced_redblue(256, cmap_gamma));
    cb = colorbar(ax, 'eastoutside');
    cb.Ticks = [cfg.climits(1), cfg.climits(2)];
    cb.TickLabels = {'Synergistic', 'Redundant'};
    cb.Label.String = 'Co-I (bits)';
    tighten_colorbar(cb);
end

function plot_label_tile(fig, tl, tile_pos, label_text, color_value, x_shift)
    ax = shifted_tile_axes(fig, tl, tile_pos, x_shift);
    axis(ax, 'off');
    text(ax, 0.5, 0.5, label_text, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', color_value, ...
        'FontWeight', 'bold', ...
        'FontSize', 12);
end

function ax = shifted_tile_axes(fig, tl, tile_pos, x_shift)
    placeholder = nexttile(tl, tile_pos);
    drawnow;
    pos = get(placeholder, 'Position');
    delete(placeholder);

    pos(1) = pos(1) + x_shift;
    ax = axes('Parent', fig, 'Position', pos);
end

function x_shift = get_left_column_shift(cfg)
    if isfield(cfg, 'left_column_shift') && ~isempty(cfg.left_column_shift)
        x_shift = cfg.left_column_shift;
    else
        x_shift = 0.1;
    end
end

function set_standard_colorbar_ticks(cb, climits, nTicks)
    if nargin < 3 || isempty(nTicks)
        nTicks = 5;
    end
    if isempty(climits) || numel(climits) ~= 2 || any(~isfinite(climits))
        return;
    end
    cb.Limits = climits;
    cb.Ticks = linspace(climits(1), climits(2), nTicks);
end

function lift_title(ax)
    th = get(ax, 'Title');
    if isempty(th) || ~isvalid(th)
        return;
    end
    try
        set(th, 'Units', 'normalized', 'Position', [0.5 1.03 0]);
    catch
    end
end

function tighten_colorbar(cb)
    try
        cb.FontSize = 7;
        cb.Label.FontSize = 8;
    catch
    end
end

function [mi1_col, mi2_col] = get_mi_colors(label)
    label_lc = lower(label);
    orange = [0.8500, 0.3250, 0.0980];
    green = [0 0.5 0];

    if contains(label_lc, 'central-peripheral') || contains(label_lc, 'central_peripheral')
        mi1_col = orange;
        mi2_col = green;
    elseif contains(label_lc, 'central')
        mi1_col = orange;
        mi2_col = orange;
    else
        mi1_col = green;
        mi2_col = green;
    end
end

function mi_title = get_mi_title(label, panel_idx)
    label_lc = lower(label);

    if contains(label_lc, 'central-peripheral') || contains(label_lc, 'central_peripheral')
        if panel_idx == 1
            mi_title = 'Central MI';
        else
            mi_title = 'Peripheral MI';
        end
    elseif contains(label_lc, 'central')
        mi_title = 'Central MI';
    else
        mi_title = 'Peripheral MI';
    end
end


function save_plot(basefold, pname, tag, label, cond)
    % Ensures save folder exists and saves figure as .fig, .png, and .svg
    out_dir = fullfile(basefold, 'Results', pname, cond);
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    fig_path = fullfile(out_dir, sprintf('%s_%s_%s_within_%s.fig', pname, tag, label, cond));
    png_path = fullfile(out_dir, sprintf('%s_%s_%s_within_%s.png', pname, tag, label, cond));
    svg_path = fullfile(out_dir, sprintf('%s_%s_%s_within_%s.svg', pname, tag, label, cond));
    saveas(gcf, fig_path, 'fig');
    try
        exportgraphics(gcf, png_path, 'Resolution', 300);
    catch
        saveas(gcf, png_path);
    end
    % Ensure vector renderer for SVG output (text remains as vector elements)
    try
        set(gcf, 'Renderer', 'painters');
    catch
    end
    try
        exportgraphics(gcf, svg_path, 'ContentType', 'vector');
    catch
        print(gcf, svg_path, '-dsvg');
    end
    close(gcf);
end

function cmap = enhanced_redblue(n, gamma)
% Increase contrast around the center while preserving full endpoints.
if nargin < 1 || isempty(n)
    n = 256;
end
if nargin < 2 || isempty(gamma)
    gamma = 0.65;
end

base = redblue(n);
t = linspace(-1, 1, n)';
tw = sign(t) .* abs(t).^gamma;
idx = 1 + (n - 1) * (tw + 1) / 2;

cmap = zeros(n, 3);
for c = 1:3
    cmap(:, c) = interp1(1:n, base(:, c), idx, 'linear');
end
end
