%% GETS THE DATA INTO MATRICES THAT ARE CONVINIENT FOR PLOTTING
%Input:
%       Datatype: 'Marmo_EcoG', i.e., type of data (human ecog, marmo etc)
%       partname: 'Ji'
%       condition: 'XY'
%       cutoff: final electrode number that is temporal. Eg. 16 if electrodes 1-16  are temporal and 17-27 are frontal
%       times: times vector for the data you are using 

%Output
% TempFFIall: All temporal electrode concaneted into one
% TempFFi: mean of all electrodes of CoI
%

function [tempFFiall, tempFFi, tempFFm, tempFFmr, tempFFms,tempFFmi1, tempFFmi2,tempFFir,tempFFis, temp_nb, front_nb] = get_plotting_CoI(datatype, partname, condition, cutoff, times)
addpath('D:\iEEG\Data\All_data\EoI')

%% Load data and set up variables

%%Get data from individual files to one struct (CoI)
[~, CoI] = load_data_monkey_JI_XY();

% Set empty matrices / zeros
%Matrix size defined by length of trial (times x times)
tempFFiall = [];
tempFFm = zeros(length(times),length(times));
tempFFmr = zeros(length(times),length(times));
tempFFms = zeros(length(times),length(times));
tempFFmi1 = [];
tempFFmi2 = [];

%Loads elecs of I from EoI-file
EoI_file = strcat('EoI_data_',datatype,'.mat');
load(EoI_file)
all_elecs =  EoI.(partname).(condition);

%% Distinguish between temporal and frontal electrodes
e_count =  1;
electrode = 1;
temp_elecs = {};

%get temporal
while electrode < cutoff
electrode = char(all_elecs(e_count));
temp_elecs(e_count) = {electrode};
electrode = str2double(electrode(2:end));
e_count = e_count + 1;
end

%get frontal
front_elecs = all_elecs(e_count:end);

%% Get data into single plottable matrices

for temp_eleci = 1:length(temp_elecs)
    elec_name = char(strcat(temp_elecs(temp_eleci),'_',temp_elecs(temp_eleci)));
    tempFFiall = squeeze(cat(3,tempFFiall,CoI.(partname).(condition).data.(elec_name)));
    tempFFm = tempFFm + logical(CoI.(partname).(condition).sigMask.(elec_name));
    tempFFmr = tempFFmr + logical(CoI.(partname).(condition).sigMask.(elec_name)>0);
    tempFFms = tempFFms + logical(CoI.(partname).(condition).sigMask.(elec_name)<0);
    tempFFmi1 = squeeze(cat(3,tempFFmi1, CoI.(partname).(condition).MI1.(elec_name)));
    tempFFmi2 = squeeze(cat(3,tempFFmi2, CoI.(partname).(condition).MI2.(elec_name)));
end

tempFFi = squeeze(mean(tempFFiall,3));
tempFFir = tempFFi;
tempFFis = tempFFi;
tempFFir(tempFFi <= 0 ) = 0; 
tempFFis(tempFFi >= 0 ) = 0; 

temp_nb = length(temp_elecs);
front_nb = length(front_elecs);
end





