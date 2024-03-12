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

function [basefold, datatype, subject, all_con, condition, participants, EoI, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(get_elec)
%% TYPE OF DATA
basefold = 'D:\iEEG\Data\All_data\';
data_index = 1;
datatypes = {'Cat_EcoG', 'Marmo_EcoG', 'Human_EcoG', 'Human_EEG', 'Human_attention'}; %change according to your datatypes
datatype = char(datatypes(data_index));

%% PARTICIPANT
subject   = 1 ;

% CAT_ECOG
participants = {'Negri'};
%MARMOSETS
% participants = {'Ji' 'Nr'};
%HUMAN_ECOG
% participants = {'HEC010_CG' 'HEC011_RE' 'HEC013_DR' 'HRM011_MC' 'HRM014_VT'};
%HUMAN_EEG
%participants = {'P1};
% Human_attention 
% participants = {'AML' 'CSJ' 'EOW' 'FOD' 'JLL' 'MOM' 'OSN' 'PAS' 'RSM' 'TMC' 'VLB'};


%% CONDITION
pick_con = 2;
%Monkeys
% all_con = {'XX' 'XY' 'XX_BB' 'XY_BB'};  %change according to your conditions
%Cats
all_con = {'XX' 'rowing'};
%Human attention
% all_con = {'attention' 'counting'};
condition = char(all_con(pick_con));

%% ELECTRODES OF INTEREST
%Import the electrodes of interest from the mat-file (Run MI_ERP to get this)
if get_elec == 1
    [EoI] = EoI_data(basefold,datatype,all_con, subject, condition);
    EoI = EoI.(char(participants(subject))).(char(condition));
else
    EoI = [];
end

%% EPOCHING
%Re-epoch or not, and the epoching parameters
re_epoch = 0; % 1 = re-epoch
if re_epoch == 1
    dev_epochs = {'200'};
    std_epochs = {'100'};
    epoch_length = [-0.2 0.4];
else
    dev_epochs = [];
    std_epochs = [];
    epoch_length = [];
end

%% FILTERING & DATA
srate = 500;%Sampling rate
low_cutoff = 1;%low-cut filter
high_cutoff = 40;%high_cut filter
% filt_order = 330;%filter order

%Define baseline normalisation, cutt-off points for trials (time period of
%interest), and number of permutations for CoI
baseline = 1:50;
start_cut_off = 1:100;
end_cut_off = 226; %:end
kperm = 10;
end
