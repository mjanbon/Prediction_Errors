function max_main(task_id)
%%MAIN Entry point of Matlab SLURM job
tic
USING_HPC = 1;
%% Step 1: define parameter settings
if USING_HPC == 1
    addpath(genpath('/home/mj649/rds/hpc-work/CNM')); % Add matlab paths to code folders and subfolders
    addpath(genpath('/home/mj649/rds/hpc-work/GCMI_master'));
    addpath(genpath('/home/mj649/rds/hpc-work/Prediction_Errors'));
end

% eeglab

%Get parameters for the analysis
[basefold, datatype, all_con, condition, subject, participants, EoI,...
    srate, activity_tag, deviant_group_number, standard_group_number, corrected, stim_onset, baseline,...
    start_cut_off, end_cut_off, kperm] = Max_get_param(USING_HPC,1);
%[basefold, datatype, subject, ~ , condition, participants, EoI, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(1);

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
[S, C, D] = ndgrid(subject, condition, permutations);
param_table = table(S(:), C(:), D(:), 'VariableNames', {'Subject', 'Condition', 'electrode_x_electrode'});
%data_folder = strcat(basefold,'Co_I');


%% Step 2: Fetch task_id from command-line
%task_id = 1;
params = table2struct(param_table(task_id, :));
%params.data_folder = data_folder;


%% Step 3: Run the actual job and save the results
% [data] = max_Get_COI(params.electrode_x_electrode,basefold, datatype,...
%     subject, all_con, condition, participants,deviant_group_number,...
%     standard_group_number,corrected, srate, baseline, kperm);
[data] = max_Get_COI(permutations(task_id),basefold, datatype,...
    subject, all_con, condition, participants,deviant_group_number,...
    standard_group_number,corrected, srate, baseline, kperm);

%participantname = participants(params.Subject);
participantname = participants(subject);
patname = char(strcat(participantname,'_', permutations(task_id)));
%results_dir = strcat(basefold,participantname,'_',condition);
results_dir = strcat('/home/mj649/rds/hpc-work/Drosophila_Results/Drosophila_CoI/',...
    participantname,'_',activity_tag, '_', condition);

results_folder =  char(results_dir);
if ~exist(results_folder, 'dir'); mkdir(results_folder); end

CoI.(patname) = data;
cd (results_folder)
save (patname, 'CoI', '-v7.3' )
toc
end