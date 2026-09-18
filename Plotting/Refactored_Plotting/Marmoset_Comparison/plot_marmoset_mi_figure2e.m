function marmosetMI = plot_marmoset_mi_figure2e(excelFile, outputDir, nOutputSamples, stimOnsetIndex)
%PLOT_MARMOSET_MI_FIGURE2E Extract and plot marmoset MI-ERP data from Figure 2e.
%
% The source workbook stores Figure 2e as two MI-ERP blocks:
%   columns 1-4: MI-ERP (Temporal), monkeys Kr, Go, Fr
%   columns 6-9: MI-ERP (Frontal),  monkeys Kr, Go, Fr
%
% Temporal is treated as the marmoset peripheral analogue.
% Frontal is treated as the marmoset central analogue.
%
% Outputs:
%   marmoset_figure2e_mi.csv
%   marmoset_figure2e_mi_timecourses.png/.fig/.svg
%
% By default, the 300 extracted samples are linearly interpolated to 450
% samples. Sample 100 is treated as stimulus onset and labelled 0 ms.

if nargin < 1 || isempty(excelFile)
    excelFile = 'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\41467_2024_48329_MOESM4_ESM.xlsx';
end

if nargin < 2 || isempty(outputDir)
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    outputDir = fullfile(repoRoot, 'Data', 'Results', 'Marmoset', 'Marmoset_MI');
end

if nargin < 3 || isempty(nOutputSamples)
    nOutputSamples = 450;
end

if nargin < 4 || isempty(stimOnsetIndex)
    stimOnsetIndex = 100;
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

sheetName = 'Figure 2e';
data = readcell_with_onedrive_fallback(excelFile, sheetName);

temporalLabelCol = find_label_column(data, 'MI-ERP (Temporal)');
frontalLabelCol = find_label_column(data, 'MI-ERP (Frontal)');

temporalCols = temporalLabelCol + (1:3);
frontalCols = frontalLabelCol + (1:3);

temporalMI = cell_block_to_numeric(data(2:end, temporalCols));
frontalMI = cell_block_to_numeric(data(2:end, frontalCols));

temporalMI = interpolate_rows(temporalMI, nOutputSamples);
frontalMI = interpolate_rows(frontalMI, nOutputSamples);

sampleIndex = (1:nOutputSamples)';
timeMs = sampleIndex - stimOnsetIndex;

temporalMean = mean(temporalMI, 2, 'omitnan');
frontalMean = mean(frontalMI, 2, 'omitnan');

marmosetMI = table( ...
    sampleIndex, timeMs, ...
    temporalMI(:,1), temporalMI(:,2), temporalMI(:,3), temporalMean, ...
    frontalMI(:,1), frontalMI(:,2), frontalMI(:,3), frontalMean, ...
    'VariableNames', {'SampleIndex', 'TimeMs', ...
    'Temporal_Kr', 'Temporal_Go', 'Temporal_Fr', 'TemporalMean_Peripheral', ...
    'Frontal_Kr', 'Frontal_Go', 'Frontal_Fr', 'FrontalMean_Central'});

csvFile = fullfile(outputDir, 'marmoset_figure2e_mi.csv');
writetable(marmosetMI, csvFile);

fig = figure('Color', 'w', 'Position', [100 100 950 420]);
hold on

plot(timeMs, temporalMI(:,1), 'Color', [0.55 0.70 0.95], 'LineWidth', 0.8, 'DisplayName', 'Temporal monkeys');
plot(timeMs, temporalMI(:,2:3), 'Color', [0.55 0.70 0.95], 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(timeMs, frontalMI(:,1), 'Color', [0.95 0.65 0.55], 'LineWidth', 0.8, 'DisplayName', 'Frontal monkeys');
plot(timeMs, frontalMI(:,2:3), 'Color', [0.95 0.65 0.55], 'LineWidth', 0.8, 'HandleVisibility', 'off');
plot(timeMs, temporalMean, 'Color', [0.10 0.35 0.85], 'LineWidth', 2.2, 'DisplayName', 'Temporal mean (Peripheral)');
plot(timeMs, frontalMean, 'Color', [0.85 0.25 0.10], 'LineWidth', 2.2, 'DisplayName', 'Frontal mean (Central)');
xline(0, 'k--', 'Stim onset', 'LabelVerticalAlignment', 'bottom');

xlabel('Time from stimulus onset (ms)');
ylabel('Mutual information');
title('Marmoset Figure 2e MI-ERP');
subtitle(sprintf('300 source samples interpolated to %d samples; sample %d labelled 0 ms', ...
    nOutputSamples, stimOnsetIndex));
legend('Location', 'best', 'Box', 'off');
grid on
box off
hold off

pngFile = fullfile(outputDir, 'marmoset_figure2e_mi_timecourses.png');
figFile = fullfile(outputDir, 'marmoset_figure2e_mi_timecourses.fig');
svgFile = fullfile(outputDir, 'marmoset_figure2e_mi_timecourses.svg');
saveas(fig, pngFile);
saveas(fig, figFile);
print(fig, svgFile, '-dsvg');

fprintf('Saved marmoset MI CSV to: %s\n', csvFile);
fprintf('Saved marmoset MI plot to: %s\n', pngFile);
end

function data = readcell_with_onedrive_fallback(excelFile, sheetName)
try
    data = readcell(excelFile, 'Sheet', sheetName);
catch firstError
    [~, name, ext] = fileparts(excelFile);
    localCopy = fullfile(tempdir, [name ext]);
    try
        copyfile(excelFile, localCopy, 'f');
        data = readcell(localCopy, 'Sheet', sheetName);
    catch secondError
        error('Could not read %s sheet %s directly or from temp copy.\nDirect error: %s\nTemp-copy error: %s', ...
            excelFile, sheetName, firstError.message, secondError.message);
    end
end
end

function col = find_label_column(data, label)
isMatch = false(size(data));
for r = 1:size(data, 1)
    for c = 1:size(data, 2)
        value = data{r, c};
        if ischar(value) || isstring(value)
            isMatch(r, c) = strcmpi(strtrim(char(value)), label);
        end
    end
end

[~, col] = find(isMatch, 1);
if isempty(col)
    error('Could not find label "%s" in sheet.', label);
end
end

function numericBlock = cell_block_to_numeric(cellBlock)
numericBlock = nan(size(cellBlock));
for r = 1:size(cellBlock, 1)
    for c = 1:size(cellBlock, 2)
        value = cellBlock{r, c};
        if isnumeric(value)
            numericBlock(r, c) = value;
        elseif ischar(value) || isstring(value)
            parsed = str2double(value);
            if ~isnan(parsed)
                numericBlock(r, c) = parsed;
            end
        end
    end
end
end

function interpolatedBlock = interpolate_rows(numericBlock, nOutputSamples)
sourceX = linspace(1, nOutputSamples, size(numericBlock, 1));
targetX = 1:nOutputSamples;
interpolatedBlock = nan(nOutputSamples, size(numericBlock, 2));

for c = 1:size(numericBlock, 2)
    sourceY = numericBlock(:, c);
    valid = isfinite(sourceY);
    if nnz(valid) < 2
        continue
    end
    interpolatedBlock(:, c) = interp1(sourceX(valid), sourceY(valid), targetX, 'linear', 'extrap');
end
end
