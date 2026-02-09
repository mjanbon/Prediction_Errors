%% Simple test script to verify SLURM access
% This script logs system information and timing to a file

output_dir = '/data/SBBS-PIDProject/Maxime/Logs/SLURM_Test/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Create output filename with timestamp
timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
output_file = fullfile(output_dir, sprintf('test_slurm_%s.txt', timestamp));

% Open file for writing
fid = fopen(output_file, 'w');

% Write test information
fprintf(fid, 'SLURM Test Run\n');
fprintf(fid, '==============\n');
fprintf(fid, 'Timestamp: %s\n', datetime('now'));
fprintf(fid, 'MATLAB Version: %s\n', version);
fprintf(fid, 'Computer: %s\n', getenv('HOSTNAME'));
fprintf(fid, 'User: %s\n', getenv('USER'));
fprintf(fid, 'Working Directory: %s\n', pwd);

% Get SLURM environment variables if available
slurm_job_id = getenv('SLURM_JOB_ID');
slurm_task_id = getenv('SLURM_ARRAY_TASK_ID');

fprintf(fid, '\nSLURM Environment Variables:\n');
fprintf(fid, 'SLURM_JOB_ID: %s\n', slurm_job_id);
fprintf(fid, 'SLURM_ARRAY_TASK_ID: %s\n', slurm_task_id);

% Test some basic MATLAB operations
fprintf(fid, '\nBasic MATLAB Test:\n');
test_array = rand(100, 100);
fprintf(fid, 'Created 100x100 random array\n');
fprintf(fid, 'Mean: %f\n', mean(test_array, 'all'));

fprintf(fid, '\nTest completed successfully!\n');
fclose(fid);

% Also display to console
disp('Test completed. Output written to:');
disp(output_file);
