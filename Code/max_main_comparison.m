function max_main_comparison(task_id)
%%MAIN Entry point of Matlab SLURM job
tic
% Define parameters
USING_HPC = 2;
max_number_permutations = 225;
THRESHOLD = 0.05;
%% Step 1: define parameter settings
if USING_HPC == 1
    addpath(genpath('/home/mj649/rds/hpc-work/CNM')); % Add matlab paths to code folders and subfolders
    addpath(genpath('/home/mj649/rds/hpc-work/GCMI_master'));
    addpath(genpath('/home/mj649/rds/hpc-work/Prediction_Errors'));
end

if USING_HPC == 2 % Add matlab paths to code folders and subfolders on QMUL Apocrita
    addpath(genpath('/data/SBBS-PIDProject/Maxime/CNM')); % Add matlab paths to code folders and subfolders
    addpath(genpath('/data/SBBS-PIDProject/Maxime/GCMI_master'));
    addpath(genpath('/data/SBBS-PIDProject/Maxime/Prediction_Errors'));
end


%Get parameters for the analysis
[basefold, datatype, all_con, condition,subject,participants, EoI,...
    srate, activity_tag_1, deviant_group_number_1, standard_group_number_1,...
    activity_tag_2, deviant_group_number_2, standard_group_number_2, corrected,...
    stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_comparison_param(USING_HPC, 0)

[S, C] = ndgrid(participants, all_con);
param_table = table(S(:), C(:), 'VariableNames', {'Subject', 'Condition'}); % Get a table of all subject condition combinations
params_subjcon = table2struct(param_table(floor((task_id-1)/max_number_permutations)+1, :));
strcat(basefold, datatype, '/', char(params_subjcon.Subject), '_', char(params_subjcon.Condition), '.mat')
if exist(strcat(basefold, datatype, '/', char(params_subjcon.Subject), '_', char(params_subjcon.Condition), '.mat'), "file") == 0
    return
end
%Load EoI (now the EoI for the comparison of two activities)
EOI_filename = strcat(basefold, 'DataEoI/', 'EoI_data','_',datatype,'.mat');
load(EOI_filename,'EoI');
EoI = EoI.(char(params_subjcon.Subject)).(strcat(char(activity_tag_1),'_',char(activity_tag_2))).(char(params_subjcon.Condition));


%Get permutations
elecs = {};
permutations = {};
for electrodi = 1:length(EoI)
    elec = EoI(electrodi);
    for perm_elec = 1:length(EoI)
        comb_perms = strcat(elec,'_permuted_',EoI(perm_elec));
        elecs(perm_elec,1) = comb_perms;      
    end
    permutations = vertcat(permutations,elecs);
end

%Make table of all the permutations
[D]=ndgrid(permutations);
%[S, C, D] = ndgrid(subject, condition, permutations);
param_table = table(D(:), 'VariableNames', {'electrode_x_electrode'});
%param_table = table(S(:), C(:), D(:), 'VariableNames', {'Subject', 'Condition', 'electrode_x_electrode'});
%data_folder = strcat(basefold,'Co_I');


%% Step 2: Fetch task_id from command-line
modded_task_id = mod(task_id-1, max_number_permutations) + 1;% Find modded task_id, giving correct line in permutation table (1-indexed)
if modded_task_id > height(param_table) % Check the row of the param table exists
    return
end
params = table2struct(param_table(modded_task_id, :));


%% Step 3: Run the actual job and save the results

[data] = max_Get_COI_comparison(params.electrode_x_electrode,basefold, datatype,...
    char(params_subjcon.Subject), char(params_subjcon.Condition),deviant_group_number_1,...
    standard_group_number_1, deviant_group_number_2, standard_group_number_2,...
    corrected, srate, baseline, kperm, THRESHOLD);

%participantname = participants(params.Subject);
participantname = char(params_subjcon.Subject);
patname = char(strcat(participantname,'_', params.electrode_x_electrode));
%results_dir = strcat(basefold,participantname,'_',condition);
%Change save name to comprise the two activity tags
if USING_HPC == 1
    results_dir = strcat('/home/mj649/rds/hpc-work/Drosophila_Results/Drosophila_CoI/Comparisons/',...
    participantname,'_',activity_tag_1, '_',activity_tag_2, '_',...
    char(params_subjcon.Condition));
elseif USING_HPC == 2
    results_dir = strcat('/data/SBBS-PIDProject/Maxime/Drosophila_Results/Drosophila_CoI/Comparisons/',...
    participantname,'_', activity_tag_1, '_', activity_tag_2activity_tag,...
    '_', char(params_subjcon.Condition));
elseif USING_HPC == 0
    results_dir = strcat(basefold, '/Drosophila_CoI/','/Comparisons/', participantname, ...
        '_',activity_tag_1, '_',activity_tag_2, '_', char(params_subjcon.Condition))
end


results_folder =  char(results_dir);
if ~exist(results_folder, 'dir'); mkdir(results_folder); end

CoI.(patname) = data;
cd (results_folder)
save (patname, 'CoI', '-v7.3' )
toc
end