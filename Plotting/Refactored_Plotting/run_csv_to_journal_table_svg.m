%RUN_CSV_TO_JOURNAL_TABLE_SVG Example runner for the CSV-to-SVG table exporter.
%
% Edit the parameters below, then run this script in MATLAB. If csv_path is
% left empty, csv_to_journal_table_svg will automatically use the newest
% CSV file in the script folder.

csv_path = 'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline\Data\Results\Group_Analysis\Post_Point_Significance\postpoint_significance_all_conditions.csv';
out_dir = "C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Paper Figures";

title_text = 'Journal Table';
precision = 3;
include_png = true;
% Example: exclude columns by name or index, e.g. {'UnwantedCol'} or [2,4]
exclude_cols = [1,2,4,5,6];

out = csv_to_journal_table_svg(csv_path, out_dir, ...
    'Title', title_text, ...
    'Precision', precision, ...
    'IncludePNG', include_png, ...
    'ExcludeColumns', exclude_cols);

disp(out);