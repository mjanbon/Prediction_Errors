%% PARAMETERS FOR COI-PIPELINE
% IMPORTS SPECIFIED PARAMETERS FOR THE PIPELINE

%This is the only file that needs to be edited to run the pipeline for
%different condition and participants

%Data has to be organised in the same way as in the example folders
% and named according to the format:
%Ji_XX_dev.set

%Ji = participant name
%XX = condition
%dev = deviant/standard

%The folder structure for data is:
%C:/.../basefold/datatype/participant/condition/participant_data.set
%E.g.
%D:/.../All_data/Marmo_EcoG/Ji/XX/Ji_XX_dev.set

%get_elec = 1, Get electrodes of interest
%get_elec = 0, do not get electrodes of interest
%% 

%function [basefold, datatype, subject, all_con, condition, participants, EoI, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(get_elec)
function [basefold, datatype,all_con, condition,subject,participants, EoI,...
    srate, deviant_group_number, standard_group_number, corrected, stim_onset, baseline,...
    start_cut_off, end_cut_off, kperm] = Max_get_param(USING_HPC, get_elec)
%% TYPE OF DATA
if USING_HPC == 1
    basefold = '/home/mj649/rds/hpc-work/Drosophila_Data/';
elseif USING_HPC == 0
    current_folder = pwd;
    basefold = strcat(current_folder, '\Data\');
end
data_index = 2;
datatypes = {'Marmo_EcoG', 'Drosophila_LFP'}; %change according to your datatypes
datatype = char(datatypes(data_index));

%% PARTICIPANT
subject   = 1 ;


%DROSOPHILA
participants = {'R230720'};

%MARMOSETS
% participants = {'Ji' 'Nr'};



%% CONDITION
pick_block = 1;
%Monkeys
% all_con = {'XX' 'XY' 'XX_BB' 'XY_BB'};  %change according to your conditions

%Drosophila
all_con = {'B2'};

condition = char(all_con(pick_block));

%% ELECTRODES OF INTEREST
%Import the electrodes of interest from the mat-file (Run MI_ERP to get this)
% if get_elec == 1
%     [EoI] = EoI_data(basefold,datatype,all_con, subject, condition);
%     EoI = EoI.(char(participants(subject))).(char(condition));
% else
%     EoI = [];
% end
if get_elec == 1
    EOI_filename = strcat(basefold, 'DataEOI\', 'EoI_data','_',datatype);
    load(EOI_filename);
    EoI = EoI.(char(participants(subject))).(char(condition));
else
    EoI = [];
end


%% EPOCHING
%Re-epoch or not, and the epoching parameters
% re_epoch = 0; % 1 = re-epoch
% if re_epoch == 1
%     dev_epochs = {'200'};
%     std_epochs = {'100'};
%     epoch_length = [-0.2 0.4];
% else
%     dev_epochs = [];
%     std_epochs = [];
%     epoch_length = [];
% end

%% FILTERING & DATA
srate = 1000; %Sampling rate
deviant_group_number = [3, 4]; %Deviant group in groupHyper for awake
% deviant_group_number = [5,11]; %Deviant group in groupHyper for sleep/wake
standard_group_number =[1, 2]; %Carrier group in groupHyper for awake
% standard_group_number =[2, 10]; %Carrier group in groupHyper for sleep/wake

corrected = 0; % Using filtered and normalised data or not
stim_onset = 25;
%low_cutoff = 1;%low-cut filter
%high_cutoff = 40;%high_cut filter
% filt_order = 330;%filter order

%Define baseline normalisation, cutt-off points for trials (time period of
%interest), and number of permutations for CoI
baseline = 1:5;
start_cut_off = 1;
end_cut_off = 95; %:end
kperm = 500;
end
