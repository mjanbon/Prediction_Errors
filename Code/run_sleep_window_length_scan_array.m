function run_sleep_window_length_scan_array(task_id)
% Wrapper to run sleep window-length scan for a single subject selected
% by SLURM_ARRAY_TASK_ID. Designed to be called from sbatch on Apocrita.
try
    if nargin < 1 || isempty(task_id)
        error('Missing task_id argument');
    end

    USING_HPC = 2; % 2 => Apocrita paths

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(repoRoot));

    if USING_HPC == 2
        addpath(genpath('/data/SBBS-PIDProject/Maxime/CNM'));
        addpath(genpath('/data/SBBS-PIDProject/Maxime/GCMI_master'));
        addpath(genpath('/data/SBBS-PIDProject/Maxime/Prediction_Errors'));
    end

    rng(0, 'twister');

    % Get participants list from parameter helper
    [~, ~, ~, ~, ~, participants] = Max_get_comparison_param(USING_HPC, 0);

    nParticipants = numel(participants);
    if task_id < 1 || task_id > nParticipants
        fprintf('task_id %d out of range (1..%d). Skipping.\n', task_id, nParticipants);
        return
    end

    % Support cell array or string array
    if iscell(participants)
        subj = participants{task_id};
    else
        subj = participants(task_id);
    end

    subj = char(subj);
    condition_tag = 'BSLEEP';

    fprintf('Running sleep window-length scan for subject %s (task %d)\n', subj, task_id);
    run_R060721_sleep_window_length_scan(subj, condition_tag);

    fprintf('Completed subject %s\n', subj);
catch ME
    fprintf('Wrapper failed: %s\n', ME.message);
    rethrow(ME);
end
end
