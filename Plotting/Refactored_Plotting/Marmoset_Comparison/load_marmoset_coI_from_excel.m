function CoI_array = load_marmoset_coI_from_excel(excel_file, sheet_name)
%LOAD_MARMOSET_COI_FROM_EXCEL Load marmoset CoI matrix from Excel file
%
% This function reads a marmoset CoI matrix from a specified Excel sheet
% and returns it as a MATLAB array.
%
% INPUTS:
%   excel_file  - Path to the Excel file containing marmoset data
%   sheet_name  - Sheet name (e.g., 'Sheet1', 'marmoset_CoI', etc.)
%
% OUTPUT:
%   CoI_array   - 2D double array containing the CoI matrix [n_regions x n_regions]
%
% USAGE:
%   CoI_array = load_marmoset_coI_from_excel('Marmoset_CoI_figures_data.xlsx', 'Sheet1');

    if nargin < 2 || isempty(sheet_name)
        sheet_name = 1;  % Default to first sheet if not specified
    end
    
    if nargin < 1 || isempty(excel_file)
        error('excel_file path is required');
    end
    
    % Check if file exists
    if ~exist(excel_file, 'file')
        error('Excel file not found: %s', excel_file);
    end
    
    % Read the data from Excel without relying on Excel COM activation.
    % Prefer non-Excel engine paths to avoid worksheet activation failures.
    try
        % Try readmatrix first (often works where readtable does not)
        CoI_array = readmatrix(excel_file, 'Sheet', sheet_name, 'UseExcel', false);
        if isempty(CoI_array) || ~isnumeric(CoI_array)
            error('Empty or non-numeric matrix from readmatrix.');
        end
    catch
        % Continue to readtable/readcell fallbacks
    end

    if ~exist('CoI_array', 'var') || isempty(CoI_array)
        try
            % Try readtable with non-Excel engine
            opts = detectImportOptions(excel_file, 'Sheet', sheet_name);
            data_table = readtable(excel_file, opts, 'Sheet', sheet_name, 'UseExcel', false);
            data_array = table2array(data_table);
            if ~isnumeric(data_array)
                error('readtable returned non-numeric data.');
            end
            valid_cols = any(~isnan(data_array), 1);
            CoI_array = data_array(:, valid_cols);
        catch
            % Final fallback: readcell + conversion to numeric matrix
            [~, ~, file_ext] = fileparts(excel_file);
            fprintf('Excel debug: file=%s | ext=%s\n', excel_file, file_ext);
            sheets_dbg = [];
            try
                [status_txt, sheets_dbg] = xlsfinfo(excel_file);
                fprintf('Excel debug: xlsfinfo status=%s | sheet_count=%d\n', char(string(status_txt)), numel(sheets_dbg));
            catch ME_dbg
                fprintf('Excel debug: xlsfinfo failed: %s\n', ME_dbg.message);
            end

            % Try to obtain sheet names from workbook XML (does not require Excel activation)
            sheets_xml = get_sheet_names_from_xlsx_xml(excel_file);
            if ~isempty(sheets_xml)
                sheets_dbg = sheets_xml;
            end

            % Resolve requested sheet to an index if possible.
            resolved_sheet = sheet_name;
            if ischar(sheet_name) || (isstring(sheet_name) && isscalar(sheet_name))
                requested = char(string(sheet_name));
                idx = [];
                if isstring(sheets_dbg)
                    sheets_dbg = cellstr(sheets_dbg(:));
                end
                if iscell(sheets_dbg)
                    idx = find(strcmpi(sheets_dbg, requested), 1);
                end
                if isempty(idx) && iscell(sheets_dbg)
                    req_norm = regexprep(requested, '[\s_\-\.]', '');
                    dbg_norm = cellfun(@(s) regexprep(char(string(s)), '[\s_\-\.]', ''), sheets_dbg, 'UniformOutput', false);
                    idx = find(strcmpi(dbg_norm, req_norm), 1);
                end
                if ~isempty(idx)
                    resolved_sheet = idx;
                end
            end

            try
                cell_data = readcell(excel_file, 'Sheet', resolved_sheet, 'UseExcel', false);
                num_data = cell_to_numeric_matrix(cell_data);
                if isempty(num_data)
                    error('No numeric data found with readcell fallback.');
                end
                CoI_array = num_data;
            catch
                available_sheets_str = get_available_sheet_names(excel_file);
                error('Failed to read Excel file. Requested sheet: %s. Available sheets: %s', char(string(sheet_name)), available_sheets_str);
            end
        end
    end
    
    % Ensure output is 2D numeric array
    if ~isnumeric(CoI_array)
        error('Loaded data is not numeric. Check Excel file format.');
    end
    
    if ndims(CoI_array) ~= 2
        error('Expected 2D array from Excel file, got %d dimensions', ndims(CoI_array));
    end
    
    fprintf('Loaded marmoset CoI matrix: %d x %d\n', size(CoI_array, 1), size(CoI_array, 2));
end

function available_sheets_str = get_available_sheet_names(excel_file)
% Try multiple APIs to list worksheet names across MATLAB versions/platforms.
    available_sheets_str = '(unable to list sheet names)';

    % Newer MATLAB API
    try
        sheets = sheetnames(excel_file);
        if ~isempty(sheets)
            available_sheets_str = sheet_list_to_str(sheets);
            return;
        end
    catch
        % Continue to fallback
    end

    % Legacy API
    try
        [~, sheets] = xlsfinfo(excel_file);
        if ~isempty(sheets) && ~(numel(sheets) == 1 && contains(lower(char(string(sheets{1}))), 'unreadable excel file'))
            available_sheets_str = sheet_list_to_str(sheets);
            return;
        end
    catch
        % Keep default if unavailable
    end

    % XML fallback for .xlsx files (no Excel activation required)
    try
        sheets_xml = get_sheet_names_from_xlsx_xml(excel_file);
        if ~isempty(sheets_xml)
            available_sheets_str = sheet_list_to_str(sheets_xml);
            return;
        end
    catch
        % Keep default
    end
end

function sheets = get_sheet_names_from_xlsx_xml(excel_file)
% Extract worksheet names by reading xl/workbook.xml from the xlsx archive.
    sheets = {};
    [~, ~, ext] = fileparts(excel_file);
    if ~strcmpi(ext, '.xlsx')
        return;
    end

    temp_dir = tempname;
    mkdir(temp_dir);
    cleanup_obj = onCleanup(@() safe_rmdir(temp_dir));

    try
        unzip(excel_file, temp_dir);
        workbook_xml = fullfile(temp_dir, 'xl', 'workbook.xml');
        if ~exist(workbook_xml, 'file')
            return;
        end
        xml_txt = fileread(workbook_xml);
        tokens = regexp(xml_txt, 'name="([^"]+)"', 'tokens');
        if ~isempty(tokens)
            sheets = cellfun(@(t) t{1}, tokens, 'UniformOutput', false);
        end
    catch
        sheets = {};
    end

    clear cleanup_obj;
end

function safe_rmdir(path_to_remove)
    if exist(path_to_remove, 'dir')
        try
            rmdir(path_to_remove, 's');
        catch
            % Ignore cleanup errors
        end
    end
end

function out = sheet_list_to_str(sheets)
% Convert sheet name collections to a printable comma-separated string.
    if isempty(sheets)
        out = '(no sheet names found)';
        return;
    end

    if ischar(sheets)
        out = sheets;
        return;
    end

    if isstring(sheets)
        sheets = cellstr(sheets(:));
    end

    if iscell(sheets)
        parts = cell(1, numel(sheets));
        for i = 1:numel(sheets)
            parts{i} = char(string(sheets{i}));
        end
        out = parts{1};
        for i = 2:numel(parts)
            out = [out ', ' parts{i}];
        end
        return;
    end

    out = char(string(sheets));
end

function num_mat = cell_to_numeric_matrix(cell_data)
% Convert mixed readcell output to a dense numeric matrix.
    num_mat = [];
    if isempty(cell_data)
        return;
    end

    n_rows = size(cell_data, 1);
    n_cols = size(cell_data, 2);
    tmp = nan(n_rows, n_cols);

    for r = 1:n_rows
        for c = 1:n_cols
            v = cell_data{r, c};
            if isnumeric(v) && isscalar(v)
                tmp(r, c) = v;
            elseif islogical(v) && isscalar(v)
                tmp(r, c) = double(v);
            elseif ischar(v) || (isstring(v) && isscalar(v))
                vv = str2double(char(string(v)));
                if ~isnan(vv)
                    tmp(r, c) = vv;
                end
            end
        end
    end

    valid_rows = any(~isnan(tmp), 2);
    valid_cols = any(~isnan(tmp), 1);
    if any(valid_rows) && any(valid_cols)
        num_mat = tmp(valid_rows, valid_cols);
    end
end
