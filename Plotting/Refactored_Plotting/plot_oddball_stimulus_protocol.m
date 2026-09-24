function plot_oddball_stimulus_protocol()
% PLOT_ODDBALL_STIMULUS_PROTOCOL
% Main-figure schematic of the Drosophila colour oddball protocol:
% four standard flashes followed by one deviant flash, with 50 ms on and
% 50 ms off between flashes. A faded continuation indicates repetition.

%% Settings
on_ms = 50;
off_ms = 50;
n_standard = 4;
n_flash = n_standard + 1;
n_repeat_preview = 2;
show_final_off = true;

green = [0.05 0.66 0.36];
blue = [0.05 0.32 0.92];
ink = [0.07 0.08 0.09];
muted = [0.34 0.36 0.39];
rule = [0.62 0.64 0.67];
track = [0.965 0.968 0.972];
track_edge = [0.84 0.85 0.87];

basefold = 'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline\Data\';
out_dir = fullfile(basefold, 'Results', 'Stimulus_Protocol');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% Derived geometry
period_ms = on_ms + off_ms;
main_total_ms = n_flash * on_ms + (n_flash - 1 + double(show_final_off)) * off_ms;
total_ms = main_total_ms + n_repeat_preview * period_ms;

flash_onsets = (0:(n_flash + n_repeat_preview - 1)) * period_ms;
row_y = [1.42 0.58];
bar_h = 0.32;
track_h = 0.48;
rounding = 0.10;

conditions = { ...
    struct('label', 'Green standards, blue deviant', 'standard_color', green, 'deviant_color', blue), ...
    struct('label', 'Blue standards, green deviant', 'standard_color', blue, 'deviant_color', green)};

%% Figure
fig = figure('Color', 'w', 'Position', [100 100 1120 450]);
ax = axes(fig);
hold(ax, 'on');
axis(ax, 'off');

% Header.
text(ax, 0, 2.00, 'Colour oddball flash sequence', ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
    'Color', ink, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(ax, total_ms, 2.00, sprintf('%d ms on  /  %d ms off', on_ms, off_ms), ...
    'FontName', 'Arial', 'FontSize', 12.5, 'FontWeight', 'bold', 'Color', ink, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

for r = 1:numel(conditions)
    draw_condition(ax, conditions{r}, row_y(r), flash_onsets, on_ms, off_ms, ...
        n_standard, n_flash, n_repeat_preview, main_total_ms, total_ms, ...
        bar_h, track_h, rounding, track, track_edge, ink, muted);
end

draw_time_axis(ax, total_ms, rule, muted, ink);
draw_duration_guides(ax, on_ms, off_ms, rule, muted);

xlim(ax, [-172 total_ms + 58]);
ylim(ax, [0 2.17]);

set(fig, 'Renderer', 'painters');
saveas(fig, fullfile(out_dir, 'oddball_stimulus_protocol.fig'));
print(fig, fullfile(out_dir, 'oddball_stimulus_protocol.png'), '-dpng', '-r300');
print(fig, fullfile(out_dir, 'oddball_stimulus_protocol.svg'), '-dsvg', '-painters');

fprintf('Saved oddball protocol figure to:\n%s\n', out_dir);

end

function draw_condition(ax, c, y, flash_onsets, on_ms, off_ms, n_standard, n_flash, ...
    n_repeat_preview, main_total_ms, total_ms, bar_h, track_h, rounding, track, track_edge, ink, muted)

% Soft rail behind the sequence.
rectangle(ax, 'Position', [0, y - 0.08, total_ms, track_h], ...
    'Curvature', 0.08, 'FaceColor', track, 'EdgeColor', track_edge, 'LineWidth', 0.35);

% Row label.
text(ax, -34, y + bar_h/2, c.label, ...
    'FontName', 'Arial', 'FontSize', 12, 'Color', ink, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

% Flashes: one full oddball cycle plus a faded start of the next cycle.
for k = 1:numel(flash_onsets)
    x0 = flash_onsets(k);
    is_preview = k > n_flash;

    if k <= n_standard || is_preview
        flash_color = c.standard_color;
        edge_color = 'none';
        line_width = 0.1;
        label_text = '';
    else
        flash_color = c.deviant_color;
        edge_color = ink;
        line_width = 0.85;
        label_text = 'D';
    end

    if is_preview
        flash_color = blend_with_white(flash_color, 0.48);
    end

    rectangle(ax, 'Position', [x0, y, on_ms, bar_h], ...
        'Curvature', rounding, 'FaceColor', flash_color, ...
        'EdgeColor', edge_color, 'LineWidth', line_width);

    if ~isempty(label_text)
        text(ax, x0 + on_ms/2, y + bar_h/2, label_text, ...
            'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold', ...
            'Color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

% Dots after the faded continuation indicate the sequence keeps repeating.
ellipsis_x = flash_onsets(n_flash + n_repeat_preview) + on_ms + off_ms * 0.35;
for d = 0:2
    plot(ax, ellipsis_x + d*12, y + bar_h/2, '.', 'Color', muted, 'MarkerSize', 13);
end

% Sequence brackets.
std_start = flash_onsets(1);
std_end = flash_onsets(n_standard) + on_ms;
dev_start = flash_onsets(n_standard + 1);
dev_end = dev_start + on_ms;
preview_start = flash_onsets(n_flash + 1);
preview_end = ellipsis_x + 2*12;
bracket_y = y + bar_h + 0.18;
draw_bracket(ax, std_start, std_end, bracket_y, '4 standards', muted, 10);
draw_bracket(ax, dev_start, dev_end, bracket_y, 'deviant', muted, 10);
draw_bracket(ax, preview_start, preview_end, bracket_y, 'repeat', muted, 10);

% Small off-period marks for the full cycle.
for k = 1:(n_flash - 1)
    off_start = flash_onsets(k) + on_ms;
    x_mid = off_start + off_ms/2;
    plot(ax, [x_mid x_mid], [y - 0.045 y - 0.005], '-', ...
        'Color', [0.72 0.73 0.75], 'LineWidth', 0.35);
end
end

function draw_bracket(ax, x1, x2, y, label_text, color, font_size)
plot(ax, [x1 x2], [y y], '-', 'Color', color, 'LineWidth', 0.5);
plot(ax, [x1 x1], [y y - 0.045], '-', 'Color', color, 'LineWidth', 0.5);
plot(ax, [x2 x2], [y y - 0.045], '-', 'Color', color, 'LineWidth', 0.5);
text(ax, (x1 + x2)/2, y + 0.050, label_text, ...
    'FontName', 'Arial', 'FontSize', font_size, 'Color', color, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

function draw_time_axis(ax, total_ms, rule, muted, ink)
axis_y = 0.21;
plot(ax, [0 total_ms], [axis_y axis_y], '-', 'Color', rule, 'LineWidth', 0.55);

tick_vals = 0:100:total_ms;
for t = tick_vals
    plot(ax, [t t], [axis_y - 0.025 axis_y + 0.025], '-', ...
        'Color', rule, 'LineWidth', 0.45);
    text(ax, t, axis_y - 0.080, sprintf('%d', t), ...
        'FontName', 'Arial', 'FontSize', 10, 'Color', muted, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

text(ax, total_ms/2, axis_y - 0.20, 'Time (ms)', ...
    'FontName', 'Arial', 'FontSize', 11, 'Color', ink, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

function draw_duration_guides(ax, on_ms, off_ms, rule, muted)
guide_y = 0.075;
draw_interval(ax, 0, on_ms, guide_y, 'flash on', rule, muted);
draw_interval(ax, on_ms, on_ms + off_ms, guide_y, 'off', rule, muted);
end

function draw_interval(ax, x1, x2, y, label_text, rule, muted)
plot(ax, [x1 x2], [y y], '-', 'Color', rule, 'LineWidth', 0.45);
plot(ax, [x1 x1], [y - 0.018 y + 0.018], '-', 'Color', rule, 'LineWidth', 0.45);
plot(ax, [x2 x2], [y - 0.018 y + 0.018], '-', 'Color', rule, 'LineWidth', 0.45);
text(ax, (x1 + x2)/2, y - 0.038, label_text, ...
    'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'Color', muted, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

function color_out = blend_with_white(color_in, amount)
color_out = color_in .* (1 - amount) + [1 1 1] .* amount;
end
