function main(task_id)
%%MAIN Entry point of Matlab SLURM job

%% Step 1: define parameter settings
addpath(genpath('/home/jma201/coi-ieeg'));
addpath(genpath('/home/jma201/iEEG'));
eeglab

%Get parameters for the analysis
[basefold, datatype, subject, ~ , condition, participants, EoI, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(1);

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
param_table = table(S(:), C(:),D(:), 'VariableNames', {'Subject', 'Condition', 'electrode_x_electrode'});
data_folder = '/home/jma201/coi-ieeg';


%% Step 2: Fetch task_id from command-line
params = table2struct(param_table(task_id, :));
params.data_folder = data_folder;


%% Step 3: Run the actual job and save the results
[data] = Get_COI(params.electrode_x_electrode,basefold, datatype, subject, condition, participants, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm);

participantname = participants(params.Subject);
patname = char(strcat(participantname,'_', params.electrode_x_electrode));
results_dir = strcat('/rds/user/jma201/hpc-work/results_',participantname,'_',condition);

results_folder =  char(results_dir);
if ~exist(results_folder, 'dir'); mkdir(results_folder); end

CoI.(patname) = data;
cd (results_folder)
save (patname, 'CoI', '-v7.3' )
end