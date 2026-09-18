function out = csv_to_journal_table_svg(csv_path, out_dir, varargin)
%CSV_TO_JOURNAL_TABLE_SVG Render a CSV file as a journal-style SVG table.
%
% Usage:
%   out = csv_to_journal_table_svg(csv_path)
%   out = csv_to_journal_table_svg(csv_path, out_dir)
%   out = csv_to_journal_table_svg(csv_path, out_dir, 'Name', Value, ...)
%
% Optional name-value pairs:
%   'Title'        - Title shown above the table (default: no title)
%   'Precision'    - Decimal places for numeric values (default: 3)
%   'IncludePNG'   - Also export a PNG preview (default: true)
%   'FontName'     - Table font (default: Helvetica)
%   'BaseName'     - Output file stem (default: derived from CSV name)
%
% The script reads the CSV into a table, formats it into display strings,
% and exports a vector SVG suitable for a manuscript or supplement.

if nargin < 1 || isempty(csv_path)
    script_dir = fileparts(mfilename('fullpath'));
    candidates = dir(fullfile(script_dir, '*.csv'));
    if isempty(candidates)
        error('No CSV files found in %s.', script_dir);
    end
    [~, idx] = max([candidates.datenum]);
    csv_path = fullfile(script_dir, candidates(idx).name);
end

if nargin < 2 || isempty(out_dir)
    out_dir = fileparts(csv_path);
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

p = inputParser;
p.addParameter('Title', '', @(x) ischar(x) || isstring(x));
p.addParameter('Precision', 3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('IncludePNG', true, @(x) islogical(x) && isscalar(x));
p.addParameter('FontName', 'Helvetica', @(x) ischar(x) || isstring(x));
p.addParameter('BaseName', '', @(x) ischar(x) || isstring(x));
% Exclude columns by name (cellstr or string) or by numeric indices
p.addParameter('ExcludeColumns', {}, @(x) iscellstr(x) || isstring(x) || isnumeric(x));
p.parse(varargin{:});
opts = p.Results;

if ~exist(csv_path, 'file')
    error('CSV file not found: %s', csv_path);
end

T = readtable(csv_path, 'PreserveVariableNames', true);
if isempty(T)
    error('CSV is empty: %s', csv_path);
end

% Remove excluded columns if requested
if ~isempty(opts.ExcludeColumns)
    cols = opts.ExcludeColumns;
    varNames = T.Properties.VariableNames;
    if isnumeric(cols)
        cols = round(cols(:)');
        cols = cols(cols >= 1 & cols <= numel(varNames));
        if ~isempty(cols)
            T(:, cols) = [];
        end
    else
        cols = cellstr(string(cols));
        toRemove = ismember(varNames, cols);
        if any(toRemove)
            T(:, toRemove) = [];
        end
    end
end

displayData = table_to_display_cells(T, opts.Precision);

if isempty(opts.BaseName)
    [~, stem] = fileparts(csv_path);
else
    stem = char(opts.BaseName);
end

svg_path = fullfile(out_dir, [stem, '.svg']);
png_path = fullfile(out_dir, [stem, '.png']);

% Compute a sensible figure size so table rows have enough vertical space
nRows = size(displayData, 1);
nCols = size(displayData, 2);
fig_w = max(800, min(1800, 160 * nCols));
fig_h = max(240, 48 * nRows + 80);
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, fig_w, fig_h]);
render_table_figure(fig, displayData, string(opts.Title), char(opts.FontName));

try
    exportgraphics(fig, svg_path, 'ContentType', 'vector');
catch
    print(fig, svg_path, '-dsvg');
end

if opts.IncludePNG
    try
        exportgraphics(fig, png_path, 'Resolution', 300);
    catch
        frame = getframe(fig);
        imwrite(frame.cdata, png_path);
    end
end

close(fig);

out = struct();
out.csv = csv_path;
out.svg = svg_path;
out.png = png_path;
out.table = T;

fprintf('Saved SVG table: %s\n', svg_path);
if opts.IncludePNG
    fprintf('Saved PNG preview: %s\n', png_path);
end

end

function displayData = table_to_display_cells(T, precision)
varNames = T.Properties.VariableNames;
nRows = height(T);
nCols = width(T);
displayData = cell(nRows + 1, nCols);
displayData(1, :) = varNames;

for c = 1:nCols
    col = T.(varNames{c});
    for r = 1:nRows
        displayData{r + 1, c} = format_value(col(r), precision);
    end
end
end

function s = format_value(x, precision)
if ismissing(x)
    s = '';
elseif isnumeric(x) || islogical(x)
    if isnumeric(x) && isscalar(x) && isnan(x)
        s = 'NaN';
    elseif islogical(x)
        s = string(x);
        s = char(s);
    else
        fmt = sprintf('%%.%df', precision);
        if isscalar(x)
            s = sprintf(fmt, x);
        else
            s = strjoin(arrayfun(@(v) sprintf(fmt, v), x(:)', 'UniformOutput', false), ', ');
        end
    end
elseif isstring(x) || ischar(x) || iscategorical(x)
    s = char(string(x));
else
    try
        s = char(string(x));
    catch
        s = '<unprintable>';
    end
end
end

function render_table_figure(fig_handle, table_data, title_text, font_name)
ax = axes(fig_handle);
axis(ax, 'off');
set(ax, 'Position', [0 0 1 1]);

rows = size(table_data, 1);
cols = size(table_data, 2);

% Estimate relative column widths from text lengths, with sensible bounds.
char_counts = zeros(1, cols);
for c = 1:cols
    for r = 1:rows
        char_counts(c) = max(char_counts(c), strlength(string(table_data{r, c})));
    end
end
char_counts = double(char_counts);
char_counts(char_counts < 6) = 6;
col_widths = char_counts ./ sum(char_counts);

left = 0.03;
right = 0.97;
top = 0.95;
bottom = 0.04;
title_gap = 0.07;

if strlength(title_text) > 0
    text(ax, 0.5, 0.985, title_text, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontName', font_name, 'FontSize', 14, ...
        'FontWeight', 'bold', 'Interpreter', 'none');
    top = top - title_gap;
end

table_h = top - bottom;
row_h = table_h / rows;
table_w = right - left;
col_widths = col_widths * table_w;
col_x = left + [0, cumsum(col_widths(1:end-1))];

for r = 1:rows
    y0 = top - r * row_h;
    if r == 1
        face = [0.89 0.89 0.92];
        weight = 'bold';
    else
        if mod(r, 2) == 0
            face = [0.98 0.98 0.98];
        else
            face = [1 1 1];
        end
        weight = 'normal';
    end

    rectangle(ax, 'Position', [left, y0, table_w, row_h], ...
        'FaceColor', face, 'EdgeColor', [0.80 0.80 0.80], 'LineWidth', 0.5);

    for c = 1:cols
        if c == 1
            x = col_x(c) + 0.008;
            halign = 'left';
        else
            x = col_x(c) + col_widths(c) / 2;
            halign = 'center';
        end
        text(ax, x, y0 + row_h / 2, table_data{r, c}, ...
            'HorizontalAlignment', halign, ...
            'VerticalAlignment', 'middle', ...
            'FontName', font_name, ...
            'FontSize', 10, ...
            'FontWeight', weight, ...
            'Interpreter', 'none');
    end
end

axis(ax, [0 1 0 1]);
set(ax, 'YDir', 'normal');
end