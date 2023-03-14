%% GETS THE DATA INTO MATRICES THAT ARE CONVINIENT FOR PLOTTING

%% INPUT

%       Datatype: 'Marmo_EcoG', i.e., type of data (human ecog, marmo etc)
%       partname: 'Ji'
%       condition: 'XY'
%       cutoff: final electrode number that is temporal. Eg. 16 if
%       electrodes 1-16  are temporal and 17-27 are frontal (or any other
%       cut off-point inherent to the electrodes
%       times: times vector for the data you are using 

%% OUTPUT

%       tempFFIall: CoI of all temporal electrodes concateneted into one
%       (within)
%       tempFFi: mean of all electrodes of CoI
%       tempFFir: only redundant
%       tempFFis: only synergetic 
%       tempFFm: masks
%       tempFFmr:redundant masks
%       tempFFms: synergetic masks

%       FrontFFI ...
%       ... =frontal within

%       tempfrontFFI ...
%       ... =temporo-frontal within

%       tempintFFI ...
%       ... =temporal interactions within

%       frontintFFI ...
%       ... =temporal interactions within
%%


function [tempFFiall, tempFFi, temp_all_mask, tempFFm, tempFFmr, tempFFms,tempFFmi1, tempFFmi2,tempFFir,tempFFis, temp_nb, front_nb, tempfront_nb, tempint_nb, frontint_nb, ...
          frontFFiall, frontFFi, front_all_mask, frontFFm, frontFFmr, frontFFms, frontFFmi1, frontFFmi2, frontFFir, frontFFis,...
          tempfrontFFiall, tempfrontFFi, tempfrontFFm, tempfrontFFmr, tempfrontFFms, tempfrontFFmi1, tempfrontFFmi2, tempfrontFFir, tempfrontFFis,...
          tempintFFiall, tempintFFi, tempintFFm, tempintFFmr, tempintFFms, tempintFFmi1, tempintFFmi2, tempintFFir, tempintFFis,...
          frontintFFiall, frontintFFi, frontintFFm, frontintFFmr, frontintFFms, frontintFFmi1, frontintFFmi2, frontintFFir, frontintFFis] = get_plotting_CoI(basefold, datatype, participantnum, partname, condition, cutoff, times)

addpath('D:\iEEG\Data\All_data\EoI')

%% Load data and set up variables
cutoff = cutoff(participantnum);
%%Get data from individual files to one struct (CoI)
[~, CoI] = load_data(basefold,datatype,condition,partname);

% Set empty matrices / zeros
%Matrix size defined by length of trial (times x times)
tempFFiall = [];
tempFFm = zeros(length(times),length(times));
tempFFmr = zeros(length(times),length(times));
tempFFms = zeros(length(times),length(times));
tempFFmi1 = [];
tempFFmi2 = [];

frontFFiall = [];
frontFFm = zeros(length(times),length(times));
frontFFmr = zeros(length(times),length(times));
frontFFms = zeros(length(times),length(times));
frontFFmi1 = [];
frontFFmi2 = [];

tempfrontFFiall = [];
tempfrontFFm = zeros(length(times),length(times));
tempfrontFFmr = zeros(length(times),length(times));
tempfrontFFms = zeros(length(times),length(times));
tempfrontFFmi1 = [];
tempfrontFFmi2 = [];

tempintFFiall = [];
tempintFFm = zeros(length(times),length(times));
tempintFFmr = zeros(length(times),length(times));
tempintFFms = zeros(length(times),length(times));
tempintFFmi1 = [];
tempintFFmi2 = [];

frontintFFiall = [];
frontintFFm = zeros(length(times),length(times));
frontintFFmr = zeros(length(times),length(times));
frontintFFms = zeros(length(times),length(times));
frontintFFmi1 = [];
frontintFFmi2 = [];

temp_all_mask = []; 
front_all_mask = [];

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

%Number of electrodes
temp_nb = length(temp_elecs);
front_nb = length(front_elecs);

%% Get data into single plottable matrices

%Temporal within
for temp_eleci = 1:length(temp_elecs)
    elec_name = char(strcat(temp_elecs(temp_eleci),'_',temp_elecs(temp_eleci)));
    tempFFiall = squeeze(cat(3,tempFFiall,CoI.(partname).(condition).data.(elec_name)));
    tempFFm = tempFFm + logical(CoI.(partname).(condition).sigMask.(elec_name));
    temp_all_mask = squeeze(cat(3,temp_all_mask,CoI.(partname).(condition).sigMask.(elec_name)));
    tempFFmr = tempFFmr + logical(CoI.(partname).(condition).sigMask.(elec_name)>0);
    tempFFms = tempFFms + logical(CoI.(partname).(condition).sigMask.(elec_name)<0);
    tempFFmi1 = squeeze(cat(2,tempFFmi1, CoI.(partname).(condition).MI1.(elec_name)));
    tempFFmi2 = squeeze(cat(2,tempFFmi2, CoI.(partname).(condition).MI2.(elec_name)));
end

tempFFi = squeeze(mean(tempFFiall,3));
tempFFir = tempFFi;
tempFFis = tempFFi;
tempFFir(tempFFi <= 0 ) = 0; 
tempFFis(tempFFi >= 0 ) = 0; 

%Frontal within 
for front_eleci = 1:length(front_elecs)
    elec_name = char(strcat(front_elecs(front_eleci),'_',front_elecs(front_eleci)));
    frontFFiall = squeeze(cat(3,frontFFiall,CoI.(partname).(condition).data.(elec_name)));
    front_all_mask = squeeze(cat(3,front_all_mask,CoI.(partname).(condition).sigMask.(elec_name)));
    frontFFm = frontFFm + logical(CoI.(partname).(condition).sigMask.(elec_name));
    frontFFmr = frontFFmr + logical(CoI.(partname).(condition).sigMask.(elec_name)>0);
    frontFFms = frontFFms + logical(CoI.(partname).(condition).sigMask.(elec_name)<0);
    frontFFmi1 = squeeze(cat(2,frontFFmi1, CoI.(partname).(condition).MI1.(elec_name)));
    frontFFmi2 = squeeze(cat(2,frontFFmi2, CoI.(partname).(condition).MI2.(elec_name)));
end

frontFFi = squeeze(mean(frontFFiall,3));
frontFFir = frontFFi;
frontFFis = frontFFi;
frontFFir(frontFFi <= 0 ) = 0; 
frontFFis(frontFFi >= 0 ) = 0; 


%Temporo-frontal
nb = 0; %count number of electrode pairs
for temp_eleci = 1:length(temp_elecs)
    for front_eleci = 1:length(front_elecs)
    elec_name = char(strcat(temp_elecs(temp_eleci),'_',front_elecs(front_eleci)));
    tempfrontFFiall = squeeze(cat(3,tempfrontFFiall,CoI.(partname).(condition).data.(elec_name)));
    tempfrontFFm = tempfrontFFm + logical(CoI.(partname).(condition).sigMask.(elec_name));
    tempfrontFFmr = tempfrontFFmr + logical(CoI.(partname).(condition).sigMask.(elec_name)>0);
    tempfrontFFms = tempfrontFFms + logical(CoI.(partname).(condition).sigMask.(elec_name)<0);
    tempfrontFFmi1 = squeeze(cat(2,tempfrontFFmi1, CoI.(partname).(condition).MI1.(elec_name)));
    tempfrontFFmi2 = squeeze(cat(2,tempfrontFFmi2, CoI.(partname).(condition).MI2.(elec_name)));
    nb = nb + 1;
    end
end

tempfrontFFi = squeeze(mean(tempfrontFFiall,3));
tempfrontFFir = tempfrontFFi;
tempfrontFFis = tempfrontFFi;
tempfrontFFir(tempfrontFFi <= 0 ) = 0; 
tempfrontFFis(tempfrontFFi >= 0 ) = 0; 
tempfront_nb = nb;

%Temporo-temporal
nb = 0;
for temp_eleci = 1:length(temp_elecs)
    for temp_eleci_2 = 1:length(temp_elecs)
        
        %if elec1 == elec2, skip (as we dont want within elecs)
        E_one = char(temp_elecs(temp_eleci));
        E_two = char(temp_elecs(temp_eleci_2));
        
        if strcmp(E_one,E_two) == 1
            continue
        end

        elec_name = char(strcat(temp_elecs(temp_eleci),'_',temp_elecs(temp_eleci_2)));
        tempintFFiall = squeeze(cat(3,tempintFFiall,CoI.(partname).(condition).data.(elec_name)));
        tempintFFm = tempintFFm + logical(CoI.(partname).(condition).sigMask.(elec_name));
        tempintFFmr = tempintFFmr + logical(CoI.(partname).(condition).sigMask.(elec_name)>0);
        tempintFFms = tempintFFms + logical(CoI.(partname).(condition).sigMask.(elec_name)<0);
        tempintFFmi1 = squeeze(cat(2,tempintFFmi1, CoI.(partname).(condition).MI1.(elec_name)));
        tempintFFmi2 = squeeze(cat(2,tempintFFmi2, CoI.(partname).(condition).MI2.(elec_name)));
        nb = nb + 1;
    end
end

tempintFFi = squeeze(mean(tempintFFiall,3));
tempintFFir = tempintFFi;
tempintFFis = tempintFFi;
tempintFFir(tempintFFi <= 0 ) = 0; 
tempintFFis(tempintFFi >= 0 ) = 0; 
tempint_nb = nb;

%Fronto-frontal
nb = 0;
for front_eleci = 1:length(front_elecs)
    for front_eleci_2 = 1:length(front_elecs)

        %if elec1 == elec2, skip (as we dont want within elecs)
        E_one = char(temp_elecs(front_eleci));
        E_two = char(temp_elecs(front_eleci_2));
        
        if strcmp(E_one,E_two) == 1
            continue
        end

        elec_name = char(strcat(front_elecs(front_eleci),'_',front_elecs(front_eleci_2)));
        frontintFFiall = squeeze(cat(3,frontintFFiall,CoI.(partname).(condition).data.(elec_name)));
        frontintFFm = frontintFFm + logical(CoI.(partname).(condition).sigMask.(elec_name));
        frontintFFmr = frontintFFmr + logical(CoI.(partname).(condition).sigMask.(elec_name)>0);
        frontintFFms = frontintFFms + logical(CoI.(partname).(condition).sigMask.(elec_name)<0);
        frontintFFmi1 = squeeze(cat(2,frontintFFmi1, CoI.(partname).(condition).MI1.(elec_name)));
        frontintFFmi2 = squeeze(cat(2,frontintFFmi2, CoI.(partname).(condition).MI2.(elec_name)));
        nb = nb + 1;
    end
end

frontintFFi = squeeze(mean(frontintFFiall,3));
frontintFFir = frontintFFi;
frontintFFis = frontintFFi;
frontintFFir(frontintFFi <= 0 ) = 0; 
frontintFFis(frontintFFi >= 0 ) = 0;
frontint_nb = nb;

end





