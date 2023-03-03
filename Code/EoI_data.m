%% LOAD EOI
% Loads the electrodes of interest from the mat-file
%To get the EoI.mat run MI_ERP.m

function [EoI_all] = EoI_data(basefold,datatype, all_con, subject, condition)
addpath('/home/jma201/coi-ieeg')
addpath('/home/jma201/iEEG')

datatype = char(datatype);
EoI_path = strcat(basefold,'EoI');
participant_path = strcat(basefold,datatype);

partfilepath = strcat(EoI_path,'/');
filename = strcat(partfilepath,'EoI_data_',datatype,'.mat');
load(filename);

folders =struct2cell(dir(participant_path));
num_par = size(folders);                                                   
partnames = folders(1,3:num_par(2));

%Find electrodes of interest from the loaded data based on
%MIEoI.mat
patienti = subject;
conditioni = condition;
        sel_patient = char(partnames(patienti));
        sel_con = char(conditioni);
        patId = sel_patient;
        EoI_all.(patId).(sel_con) = EoI.(sel_patient).(sel_con);
end



