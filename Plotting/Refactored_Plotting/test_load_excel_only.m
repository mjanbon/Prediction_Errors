% test_load_excel_only.m
% Minimal script to test loading an Excel workbook.

clear; clc;

excel_file = 'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\41467_2024_48329_MOESM4_ESM.xlsx';
sheet_to_load = 1;   % Change to e.g. 'Figure 3a' if needed

if ~isfile(excel_file)
    error('Excel file not found: %s', excel_file);
end

fprintf('Excel file found:\n%s\n\n', excel_file);

% Step 1: list sheet names
try
    sheets = sheetnames(excel_file);
    fprintf('Workbook has %d sheets.\n', numel(sheets));
    disp(sheets);
catch ME
    fprintf('Could not list sheets with sheetnames: %s\n', ME.message);
    [status_txt, sheets] = xlsfinfo(excel_file);
    fprintf('xlsfinfo status: %s\n', char(string(status_txt)));
    if isempty(sheets)
        fprintf('xlsfinfo did not return sheet names.\n');
    else
        fprintf('Workbook has %d sheets (xlsfinfo).\n', numel(sheets));
        disp(sheets);
    end
end

% Step 2: load one sheet
fprintf('\nLoading sheet: %s\n', char(string(sheet_to_load)));
try
    data = readcell(excel_file, 'Sheet', sheet_to_load, 'UseExcel', false);
    fprintf('Loaded successfully with readcell (UseExcel=false).\n');
catch ME
    fprintf('readcell failed: %s\n', ME.message);
    data = readcell(excel_file, 'Sheet', sheet_to_load);
    fprintf('Loaded successfully with readcell (default engine).\n');
end

fprintf('Loaded cell array size: %d x %d\n', size(data, 1), size(data, 2));

% Keep result in workspace for inspection
assignin('base', 'excel_test_data', data);
assignin('base', 'excel_test_sheets', sheets);
fprintf('Variables exported to base workspace: excel_test_data, excel_test_sheets\n');
