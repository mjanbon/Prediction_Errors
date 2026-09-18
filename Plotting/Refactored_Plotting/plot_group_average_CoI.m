function plot_group_average_CoI(fly_list, basefold, datatype, activity_tag, condition, cutoff, times)
% PLOT_GROUP_AVERAGE_COI
% Computes and plots the average co-information (CoI) and percentage of significant timepoints
% across all flies and all regions (central, peripheral, cross-region).

regions = {'central', 'peripheral', 'central_peripheral', 'central_central', 'peripheral_peripheral'};
n_flies = numel(fly_list);
avg = struct();
resolved_time_len = [];

for i = 1:n_flies
    pname = fly_list{i};
    fprintf('Loading CoI for %s...\n', pname);
    CoI = max_get_plotting_CoI_rf(basefold, datatype, pname, activity_tag, condition, cutoff, times);

    for r = 1:numel(regions)
        region = regions{r};
        if isfield(CoI, region) && isfield(CoI.(region), 'mean') && ~isempty(CoI.(region).mean)
            s = size(CoI.(region).mean);
            disp([region ' size for ' pname ': ' mat2str(s)]);
            if numel(s) == 2 && s(1) == s(2)
                if isempty(resolved_time_len)
                    resolved_time_len = s(1);
                    for rr = 1:numel(regions)
                        reg = regions{rr};
                        avg.(reg).mean = zeros(resolved_time_len);
                        avg.(reg).mask = zeros(resolved_time_len);
                        avg.(reg).n_valid = 0;
                    end
                end
                if s(1) == resolved_time_len
                    avg.(region).mean = avg.(region).mean + CoI.(region).mean;
                    avg.(region).mask = avg.(region).mask + logical(CoI.(region).mask);
                    avg.(region).n_valid = avg.(region).n_valid + 1;
                else
                    fprintf('Skipping region %s for %s (size mismatch)\n', region, pname);
                end
            end
        end
    end
end

if isempty(resolved_time_len)
    error('No valid data found to determine time length.');
end

if nargin < 7 || isempty(times)
    warning('Time vector not provided. Defaulting to 1:N.');
    times = 1:resolved_time_len;
elseif length(times) ~= resolved_time_len
    warning('Time vector length does not match data size. Using 1:N.');
    times = 1:resolved_time_len;
end

for r = 1:numel(regions)
    region = regions{r};
    if avg.(region).n_valid > 0
        avg.(region).mean = avg.(region).mean / avg.(region).n_valid;
        avg.(region).mask_pct = 100 * avg.(region).mask / avg.(region).n_valid;
    else
        avg.(region).mean = nan(resolved_time_len);
        avg.(region).mask_pct = nan(resolved_time_len);
    end
end

fig = figure('Position', [100 100 1400 500], 'Renderer', 'painters');
tiledlayout(2,6,'TileSpacing','compact');

plot_coi_sign_key([-0.01 0.01]);

for r = 1:numel(regions)
    region = regions{r};

    ax1 = nexttile(r + 1);
    contourf(ax1, times, times, avg.(region).mean, 50, 'linecolor', 'none');
    colormap(ax1, redblue(256)); shading(ax1, 'flat'); clim(ax1, [-0.01 0.01]);
    title(ax1, [strrep(region,'_','-') ' avg CoI']); xlabel(ax1, 'Time (ms)'); ylabel(ax1, 'Time (ms)');
    daspect(ax1, [1 1 1]);
    cb1 = colorbar(ax1);
    format_coi_colorbar(cb1, [-0.01 0.01]);

    ax2 = nexttile(r + 7);
    imagesc(ax2, times, times, avg.(region).mask_pct); set(ax2,'YDir','normal');
    colormap(ax2, flipud(bone(256))); clim(ax2, [0 100]);
    cb2 = colorbar(ax2);
    cb2.Label.String = '% significant';
    title(ax2, [strrep(region,'_','-') ' % significant']); xlabel(ax2, 'Time (ms)'); ylabel(ax2, 'Time (ms)');
    daspect(ax2, [1 1 1]);
end

% Add overall title with activity tag included
sgtitle(sprintf('Group average CoI for %s - %s across %d flies', activity_tag, condition, n_flies));

% Create Fly_averages directory in Results folder if it doesn't exist
save_dir = fullfile(basefold, 'Results', 'Fly_averages');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('Created directory: %s\n', save_dir);
end

% Generate a descriptive filename
filename_base = sprintf('Group_Avg_CoI_%s_%s_%dFlies', activity_tag, condition, n_flies);

% Save as .fig (MATLAB figure format)
fig_filename = fullfile(save_dir, [filename_base '.fig']);
saveas(fig, fig_filename, 'fig');

% Save as .png (for easy viewing outside MATLAB)
png_filename = fullfile(save_dir, [filename_base '.png']);
saveas(fig, png_filename, 'png');

fprintf('Saved group average figures to:\n%s\n%s\n', fig_filename, png_filename);
end

function format_coi_colorbar(cb, climits)
    cb.Ticks = [climits(1), 0, climits(2)];
    cb.TickLabels = {sprintf('%.3g Synergistic', climits(1)), '0', sprintf('%.3g Redundant', climits(2))};
    cb.Label.String = 'Co-I (bits)';
end

function plot_coi_sign_key(climits)
    ax = nexttile(1);
    axis(ax, 'off');
    clim(ax, climits);
    colormap(ax, redblue(256));
    cb = colorbar(ax, 'eastoutside');
    cb.Ticks = [climits(1), climits(2)];
    cb.TickLabels = {'Synergistic', 'Redundant'};
    cb.Label.String = 'Co-I (bits)';
end
