function plot_marmoset_central_coi_like_flies(central_csv_file, basefold, condition, activity_tag, stim_onset_index, timing)
%PLOT_MARMOSET_CENTRAL_COI_LIKE_FLIES Plot marmoset central CoI like fly group plots.
%
% compare_drosophila_marmoset_coI loads the marmoset matrices directly from
% CSV files. This helper takes the full central-electrode CSV matrix, wraps
% it in the small region struct expected by plot_CoI_region, and saves the
% standard tiled CoI figure. Unlike compare_drosophila_marmoset_coI, this
% does not crop the marmoset matrix into the Drosophila -25:74 ms window.
%
% Example:
%   plot_marmoset_central_coi_like_flies( ...
%       'C:\path\to\marmoset_central.csv', ...
%       'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline', ...
%       [], 'marmoset', 100);
%
% Figures are saved under:
%   <basefold>/Results/Marmoset/Marmoset_CoI

    this_dir = fileparts(mfilename('fullpath'));
    addpath(this_dir);
    addpath(fileparts(this_dir));
    set_default_plotting();

    if nargin < 2 || isempty(basefold)
        basefold = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    end
    if nargin < 3 || isempty(condition)
        condition = 'Marmoset_CoI';
    end
    if nargin < 4 || isempty(activity_tag)
        activity_tag = 'central';
    end
    if nargin < 1 || isempty(central_csv_file)
        error('Provide the path to the marmoset central CSV file.');
    end
    if ~isfile(central_csv_file)
        error('Could not find central CSV file: %s', central_csv_file);
    end

    full_matrix = read_numeric_csv_matrix(central_csv_file);
    central_matrix = full_matrix;

    if size(central_matrix, 1) ~= size(central_matrix, 2)
        error('Expected a square CoI matrix, but got %d x %d.', size(central_matrix, 1), size(central_matrix, 2));
    end

    if nargin < 6 || isempty(timing)
        if nargin >= 5 && ~isempty(stim_onset_index)
            timing = (1:size(central_matrix, 1)) - round(stim_onset_index);
        else
            timing = 1:size(central_matrix, 1);
        end
    end

    if numel(timing) ~= size(central_matrix, 1)
        error('Timing vector has %d points, but the marmoset matrix is %d x %d.', ...
            numel(timing), size(central_matrix, 1), size(central_matrix, 2));
    end

    region = make_region_struct_for_plotting(central_matrix);

    tick_values = make_axis_ticks(timing);
    plot_cfg = struct( ...
        'xlimits', [min(timing) max(timing)], ...
        'ylimits_CoI', [min(timing) max(timing)], ...
        'ylimit_MI', [-0.005 0.1], ...
        'xticks_CoI', tick_values, ...
        'yticks_CoI', tick_values, ...
        'yticks_MI', 0:0.05:0.2, ...
        'x_labels', tick_values, ...
        'y_labels_CoI', tick_values, ...
        'y_labels_MI', 0:0.05:0.1, ...
        'climits', symmetric_climits(central_matrix), ...
        'colormap_gamma', 0.65, ...
        'climits_mask', [0 1]);

    plot_CoI_region(region, basefold, 'Marmoset', activity_tag, condition, timing, plot_cfg, 'Central');
end

function matrix_data = read_numeric_csv_matrix(csv_file)
    try
        matrix_data = readmatrix(csv_file);
        if isnumeric(matrix_data) && ~isempty(matrix_data)
            return;
        end
    catch ME
        fprintf('readmatrix failed for %s: %s\nTrying fopen/text parsing fallback.\n', csv_file, ME.message);
    end

    fid = fopen(csv_file, 'rt');
    if fid < 0
        error(['Could not open CSV file with readmatrix or fopen: %s\n' ...
            'If this file is in OneDrive, right-click it and choose "Always keep on this device", ' ...
            'then try again.'], csv_file);
    end

    cleanup_obj = onCleanup(@() fclose(fid));
    rows = {};
    expected_cols = [];

    while true
        this_line = fgetl(fid);
        if ~ischar(this_line)
            break;
        end
        if isempty(strtrim(this_line))
            continue;
        end

        values = str2double(strsplit(strtrim(this_line), ','));
        if isempty(expected_cols)
            expected_cols = numel(values);
        elseif numel(values) ~= expected_cols
            error('CSV row has %d columns, expected %d in file %s.', numel(values), expected_cols, csv_file);
        end

        rows{end + 1, 1} = values; %#ok<AGROW>
    end

    if isempty(rows)
        error('No numeric rows found in CSV file: %s', csv_file);
    end

    matrix_data = vertcat(rows{:});
    clear cleanup_obj;
end

function region = make_region_struct_for_plotting(mean_matrix)
    n_time = size(mean_matrix, 1);
    finite_mask = isfinite(mean_matrix);

    region = struct();
    region.mean = mean_matrix;
    region.redundant = mean_matrix .* (mean_matrix > 0);
    region.synergetic = mean_matrix .* (mean_matrix < 0);
    region.mask = double(finite_mask);
    region.mask_r = double(finite_mask & mean_matrix > 0);
    region.mask_s = double(finite_mask & mean_matrix < 0);
    region.mi1 = [];
    region.mi2 = [];
    region.FFi_all = reshape(mean_matrix, n_time, n_time, 1);
    region.n = 1;
end

function ticks = make_axis_ticks(timing)
    min_t = min(timing);
    max_t = max(timing);
    if min_t == max_t
        ticks = min_t;
        return;
    end

    raw_ticks = linspace(min_t, max_t, 5);
    ticks = unique(round(raw_ticks));
end

function climits = symmetric_climits(data)
    vals = data(isfinite(data));
    if isempty(vals)
        climits = [-0.01 0.01];
        return;
    end

    max_abs = max(abs(vals));
    if max_abs == 0
        max_abs = 0.01;
    end
    climits = [-max_abs max_abs];
end
