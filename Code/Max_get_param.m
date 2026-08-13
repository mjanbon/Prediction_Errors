function [basefold, datatype, all_con, condition, subject, participants, EoI, ...
    srate, activity_tag, deviant_group_number, standard_group_number, corrected, stim_onset, baseline, ...
    start_cut_off, end_cut_off, kperm] = Max_get_param(USING_HPC, get_elec)
% MAX_GET_PARAM   Retrieve parameters for CoI pipeline analysis
%
% SYNTAX:
%   [basefold, datatype, all_con, condition, subject, participants, EoI, ...
%    srate, activity_tag, deviant_group_number, standard_group_number, ...
%    corrected, stim_onset, baseline, start_cut_off, end_cut_off, kperm] = ...
%    Max_get_param(USING_HPC, get_elec)
%
% INPUT PARAMETERS:
%   USING_HPC: Whether Queen Mary or Cambridge HPC is used or the analysis is
%              to be run locally           
%   get_elec: - get_elec = 1: Get electrodes of interest
%             - get_elec = 0: Do not get electrodes of interest              
%
% OUTPUT PARAMETERS:
%   basefold: The base folder where the LFP data is found               
%   datatype: What animal the LFP data to be analysed corresponds to                 
%   all_con: A list of conditions the data could have been recorded under                 
%   condition: The condition to be analysed             
%   subject: Used for selecting a particular participant                 
%   participants: List of participants (different flies)           
%   EoI: List of electrode numbers with significant mutual information                      
%   srate: The sampling rate of the LFPs in Hz                   
%   activity_tag: String to show the fly activity at time of recording (eg sleep, wake etc)             
%   deviant_group_number: Rows in GroupHyper (data from Drosophila preprocessing)
%                         that correspond to deviant stimuli with the correct fly activity.     
%   standard_group_number: Rows in GroupHyper (data from Drosophila preprocessing)
%                          that correspond to standard stimuli with the correct fly activity.     
%   corrected: Whether any normalisation has been done on the ERPs in the preprocessing.                 
%   stim_onset:                 
%   baseline:                   
%   start_cut_off:              
%   end_cut_off:                
%   kperm: Number of permutations in the permutation analysis                     
%
% DESCRIPTION:
%   This is the only file that needs to be edited to run the pipeline for
%   different conditions and participants (with no condition comparisons).
  

%% TYPE OF DATA
if USING_HPC == 1
    basefold = '/home/mj649/rds/hpc-work/Drosophila_Data/';
elseif USING_HPC == 2
    basefold = '/data/SBBS-PIDProject/Maxime/Drosophila_Data/';
elseif USING_HPC == 0
    current_folder = pwd;
    %basefold = strcat(current_folder, '\Data\');
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
% participants = {'R070721','R070721'};
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
activity_tag = 'beginning_sleep';
if strcmp(activity_tag, 'sleep') == 1
    deviant_group_number = [3, 4]; %Deviant group in groupHyper for decomposed sleep/wake (all sleep)
    standard_group_number = [9, 10]; %Carrier group in groupHyper for decomposed sleep/wake (all sleep)
    srate = 1000;
elseif strcmp(activity_tag, 'mid_sleep') == 1
    deviant_group_number = [43, 44]; %Deviant group in groupHyper for decomposed sleep/wake (mid mins sleep)
    standard_group_number = [41, 42]; %Carrier group in groupHyper for decomposed sleep/wake (mid mins sleep)
    srate = 1000;
elseif strcmp(activity_tag, 'beginning_sleep') == 1
    deviant_group_number = [27, 28]; %Deviant group in groupHyper for decomposed sleep/wake (first 2 mins sleep)
    standard_group_number = [21, 22]; %Carrier group in groupHyper for decomposed sleep/wake (first 2 mins sleep)
    srate = 1000;
elseif strcmp(activity_tag, 'wake_in_sleep_dataset') == 1
    deviant_group_number = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
    standard_group_number = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)
    srate = 1000;
elseif strcmp(activity_tag, 'only_wake_dataset') == 1
    deviant_group_number = [3, 4]; %Deviant group in groupHyper for only wake dataset (wake)
    standard_group_number = [1, 2]; %Carrier group in groupHyper for only wake dataset (wake)
    srate = 1000;
elseif strcmp(activity_tag, 'only_wake_edited_v1') == 1
    deviant_group_number = [5, 6]; %Deviant group in groupHyper for only wake dataset (wake)
    standard_group_number = [1, 2]; %Carrier group in groupHyper for only wake dataset (wake)
    srate = 1000;
elseif strcmp(activity_tag, 'only_wake_oddball_11_point_5') == 1
    deviant_group_number = [3, 4]; %Deviant group in groupHyper for only wake dataset (wake)
    standard_group_number = [1, 2]; %Carrier group in groupHyper for only wake dataset (wake)
    srate = 1000;
elseif strcmp(activity_tag, 'wake') == 1
    deviant_group_number = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
    standard_group_number = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)
    srate = 1000;
elseif strcmp(activity_tag, 'wake_200') == 1
    deviant_group_number = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake at 200 Hz resampling frequency (wake)
    standard_group_number = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake at 200 Hz resampling frequency (wake)
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
