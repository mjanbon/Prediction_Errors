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
    srate, activity_tag, deviant_group_number, standard_group_number, corrected, stim_onset, baseline,...
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
participants = {'R280721'};

%MARMOSETS
% participants = {'Ji' 'Nr'};



%% CONDITION
pick_block = 1;
%Monkeys
% all_con = {'XX' 'XY' 'XX_BB' 'XY_BB'};  %change according to your conditions

%Drosophila
all_con = {'BSLEEP'};

condition = char(all_con(pick_block));

%% FILTERING & DATA
activity_tag = 'wake_200';
if strcmp(activity_tag, 'mid_sleep') == 1
    deviant_group_number = [43, 44]; %Deviant group in groupHyper for decomposed sleep/wake (mid mins sleep)
    standard_group_number = [41, 42]; %Carrier group in groupHyper for decomposed sleep/wake (mid mins sleep)
    srate = 1000;
elseif strcmp(activity_tag, 'beginning_sleep') == 1
    deviant_group_number = [27, 28]; %Deviant group in groupHyper for decomposed sleep/wake (first 2 mins sleep)
    standard_group_number = [21, 22]; %Carrier group in groupHyper for decomposed sleep/wake (first 2 mins sleep)
    srate = 1000;
elseif strcmp(activity_tag, 'only_wake_dataset') == 1
    deviant_group_number = [3, 4]; %Deviant group in groupHyper for only wake dataset (wake)
    standard_group_number = [1, 2]; %Carrier group in groupHyper for only wake dataset (wake)
    srate = 1000;
elseif strcmp(activity_tag, 'wake') == 1
    deviant_group_number = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
    standard_group_number = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)
    srate = 1000;
elseif strcmp(activity_tag, 'wake_200') == 1
    deviant_group_number = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
    standard_group_number = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)
    srate = 200;
end


% srate = 200; %Sampling rate
% %deviant_group_number = [3, 4]; %Deviant group in groupHyper for awake
% % deviant_group_number = [5,11]; %Deviant group in groupHyper for sleep/wake
% deviant_group_number = [43,44]; %Deviant group in groupHyper for decomposed sleep/wake (mid mins sleep)
% 
% %standard_group_number = [1, 2]; %Carrier group in groupHyper for awake
% % standard_group_number =[2, 10]; %Carrier group in groupHyper for sleep/wake
% standard_group_number = [41, 42]; %Carrier group in groupHyper for decomposed sleep/wake (mid mins sleep)

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

%% ELECTRODES OF INTEREST
%Import the electrodes of interest from the mat-file (Run MI_ERP to get this)
% if get_elec == 1
%     [EoI] = EoI_data(basefold,datatype,all_con, subject, condition);
%     EoI = EoI.(char(participants(subject))).(char(condition));
% else
%     EoI = [];
% end
if get_elec == 1
    EOI_filename = strcat(basefold, 'DataEoI/', 'EoI_data','_',datatype,'.mat');
    load(EOI_filename,'EoI');
    EoI = EoI.(char(participants(subject))).(char(activity_tag)).(char(condition));
else
    EoI = [];
end
end
