%% PARAMETERS FOR COI-PIPELINE
% IMPORTS SPECIFIED PARAMETERS FOR THE PIPELINE

%This is the file that needs to be edited to run the pipeline for
% comparisons across conditions (eg wake vs sleep).
% It sets two activity tags and two sets of deviant and standard group numbers.
% See Max_get_param.m for a description of the parameters.


%get_elec = 1, Get electrodes of interest
%get_elec = 0, do not get electrodes of interest
%% 

%function [basefold, datatype, subject, all_con, condition, participants, EoI, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(get_elec)
function [basefold, datatype,all_con, condition,subject,participants, EoI,...
    srate, activity_tag_1, deviant_group_number_1, standard_group_number_1,...
    activity_tag_2, deviant_group_number_2, standard_group_number_2, corrected,...
    stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_comparison_param(USING_HPC, get_elec)
%% TYPE OF DATA
if USING_HPC == 1
    basefold = '/home/mj649/rds/hpc-work/Drosophila_Data/';
elseif USING_HPC == 2
    basefold = '/data/SBBS-PIDProject/Maxime/Drosophila_Data/';
elseif USING_HPC == 0
    basefold = 'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline\Data\'
end
data_index = 2;
datatypes = {'Marmo_EcoG', 'Drosophila_LFP'}; %change according to your datatypes
datatype = char(datatypes(data_index));

%% PARTICIPANT
subject   = 1 ;


%DROSOPHILA
% participants = {'R230720','R040820_1','R040820_2', 'R050820_1'...
%     , 'R060820_1', 'R060820_2', 'R060721'};
participants = {'R060721','R070721','R080721','R150721','R210721','R220721', 'R280721', 'R290721'};
% participants = {'R060721'};
%MARMOSETS
% participants = {'Ji' 'Nr'}



%% CONDITION
pick_block = 1;
%Monkeys
% all_con = {'XX' 'XY' 'XX_BB' 'XY_BB'};  %change according to your conditions

%Drosophila
%all_con = {'B2','B2WAKE', 'B3WAKE', 'BSLEEP'};
all_con = {'BSLEEP'};

condition = char(all_con(pick_block));

%% FILTERING & DATA
% Set activity tags and deviant and standard groups in groupHyper for two conditions.
activity_tag_1 = 'wake' %First activity tag
if strcmp(activity_tag_1, 'sleep') == 1
    deviant_group_number_1 = [3, 4]; %Deviant group in groupHyper for decomposed sleep/wake (all sleep)
    standard_group_number_1 = [9, 10]; %Carrier group in groupHyper for decomposed sleep/wake (all sleep)
    srate = 1000;
elseif strcmp(activity_tag_1, 'mid_sleep') == 1
    deviant_group_number_1 = [43, 44]; %Deviant group in groupHyper for decomposed sleep/wake (mid mins sleep)
    standard_group_number_1 = [41, 42]; %Carrier group in groupHyper for decomposed sleep/wake (mid mins sleep)
    srate = 1000;
elseif strcmp(activity_tag_1, 'beginning_sleep') == 1
    deviant_group_number_1 = [27, 28]; %Deviant group in groupHyper for decomposed sleep/wake (first 2 mins sleep)
    standard_group_number_1 = [21, 22]; %Carrier group in groupHyper for decomposed sleep/wake (first 2 mins sleep)
    srate = 1000;
elseif strcmp(activity_tag_1, 'wake_in_sleep_dataset') == 1
    deviant_group_number_1 = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
    standard_group_number_1 = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)
    srate = 1000;
elseif strcmp(activity_tag_1, 'only_wake_dataset') == 1
    deviant_group_number_1 = [3, 4]; %Deviant group in groupHyper for only wake dataset (wake)
    standard_group_number_1 = [1, 2]; %Carrier group in groupHyper for only wake dataset (wake)
    srate = 1000;
elseif strcmp(activity_tag_1, 'wake') == 1
    deviant_group_number_1 = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
    standard_group_number_1 = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)
    srate = 1000;
end

activity_tag_2 = 'sleep'; %Second activity tag
if strcmp(activity_tag_2, 'sleep') == 1
    deviant_group_number_2 = [3, 4]; %Deviant group in groupHyper for decomposed sleep/wake (all sleep)
    standard_group_number_2 = [9, 10]; %Carrier group in groupHyper for decomposed sleep/wake (all sleep)
    srate = 1000;
elseif strcmp(activity_tag_2, 'mid_sleep') == 1
    deviant_group_number_2 = [43, 44]; %Deviant group in groupHyper for decomposed sleep/wake (mid mins sleep)
    standard_group_number_2 = [41, 42]; %Carrier group in groupHyper for decomposed sleep/wake (mid mins sleep)
    srate = 1000;
elseif strcmp(activity_tag_2, 'beginning_sleep') == 1
    deviant_group_number_2 = [27, 28]; %Deviant group in groupHyper for decomposed sleep/wake (first 2 mins sleep)
    standard_group_number_2 = [21, 22]; %Carrier group in groupHyper for decomposed sleep/wake (first 2 mins sleep)
    srate = 1000;
elseif strcmp(activity_tag_2, 'wake_in_sleep_dataset') == 1
    deviant_group_number_2 = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
    standard_group_number_2 = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)
    srate = 1000;
elseif strcmp(activity_tag_2, 'only_wake_dataset') == 1
    deviant_group_number_2 = [3, 4]; %Deviant group in groupHyper for only wake dataset (wake)
    standard_group_number_2 = [1, 2]; %Carrier group in groupHyper for only wake dataset (wake)
    srate = 1000;
elseif strcmp(activity_tag_2, 'wake') == 1
    deviant_group_number_2 = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
    standard_group_number_2 = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)
    srate = 1000;
end


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
