%% GETS THE DATA INTO MATRICES THAT ARE CONVINIENT FOR PLOTTING

%% INPUT

%       Datatype: 'Marmo_EcoG', i.e., type of data (human ecog, marmo etc)
%       cur_par: 'Ji'
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


function [tempFFiall, tempFFi, temp_all_mask, tempFFm, tempFFmr, tempFFms,tempFFmi1, tempFFmi2,tempFFir,tempFFis, final_temp_nb, final_front_nb, final_tempfront, final_tempint, final_frontint, ...
    frontFFiall, frontFFi, front_all_mask, frontFFm, frontFFmr, frontFFms, frontFFmi1, frontFFmi2, frontFFir, frontFFis,...
    tempfrontFFiall, tempfrontFFi, tempfrontFFm, tempfrontFFmr, tempfrontFFms, tempfrontFFmi1, tempfrontFFmi2, tempfrontFFir, tempfrontFFis,...
    tempintFFiall, tempintFFi, tempintFFm, tempintFFmr, tempintFFms, tempintFFmi1, tempintFFmi2, tempintFFir, tempintFFis,...
    frontintFFiall, frontintFFi, frontintFFm, frontintFFmr, frontintFFms, frontintFFmi1, frontintFFmi2, frontintFFir, frontintFFis, ...
    temp_elecs, front_elecs] = get_plotting_CoI(basefold, datatype, participantnum, cur_par, condition, cut, times)

[basefold, datatype, participanti, ~, ~, participants, ~, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, ~, start_cut_off, end_cut_off, kperm] = Get_param(0);
participants = {'AML' 'CSJ' 'EOW' 'FOD' 'MOM' 'PAS' 'VLB'  };
addpath('D:\iEEG\Data\All_data\EoI')

size(times)
tempFFiall = [];
tempFFm = zeros(349,349);
tempFFmr =  zeros(349,349);
tempFFms =  zeros(349,349);
tempFFmi1 = [];
tempFFmi2 = [];

frontFFiall = [];
frontFFm =  zeros(349,349);
frontFFmr =  zeros(349,349);
frontFFms =  zeros(349,349);
frontFFmi1 = [];
frontFFmi2 = [];

tempfrontFFiall = [];
tempfrontFFm =  zeros(349,349);
tempfrontFFmr =  zeros(349,349);
tempfrontFFms =  zeros(349,349);
tempfrontFFmi1 = [];
tempfrontFFmi2 = [];

tempintFFiall = [];
tempintFFm =  zeros(349,349);
tempintFFmr = zeros(349,349);
tempintFFms = zeros(349,349);
tempintFFmi1 = [];
tempintFFmi2 = [];

frontintFFiall = [];
frontintFFm =  zeros(349,349);
frontintFFmr =  zeros(349,349);
frontintFFms =  zeros(349,349);
frontintFFmi1 = [];
frontintFFmi2 = [];

temp_all_mask = [];
front_all_mask = [];


final_temp_nb = 0;
final_front_nb = 0;
final_tempint = 0;
final_frontint = 0;
final_tempfront = 0;

for participanti = 1:length(participants)

cur_par = char(participants(participanti));
%% Load data and set up va

%%Get data from individual files to one struct (CoI)
[~, CoI] = load_data(basefold,datatype,condition,cur_par);
% Set empty matrices / zeros
%Matrix size defined by length of trial (times x times)

%Loads elecs of I from EoI-file
EoI_file = strcat('EoI_data_',datatype,'.mat');
load(EoI_file)
all_elecs =  EoI.(cur_par).(condition)
cut1 = EoI.(char(participants(participanti))).cutoff_temp.attention +1
cut2 = EoI.(char(participants(participanti))).cutoff_temp.counting +1

if strcmp(condition,'attention') == 1 
cutoff = cut1;
else
cutoff = cut2;
end



%% Distinguish between temporal and frontal electrodes
e_count =  1;
electrode = 1;
temp_elecs = {};

%get temporal
while e_count < cutoff
    electrode = char(all_elecs(e_count));
    temp_elecs(e_count) = {electrode}
    electrode = str2double(electrode(2:end))
    e_count = e_count + 1
end

%get frontal
front_elecs = all_elecs(e_count:end)

%Number of electrodes
temp_nb = length(temp_elecs);
final_temp_nb = final_temp_nb + temp_nb
front_nb = length(front_elecs);
final_front_nb = final_front_nb + front_nb

%% Get data into single plottable matrices

%Temporal within
for temp_eleci = 1:length(temp_elecs)
    elec_name = char(strcat(temp_elecs(temp_eleci),'_',temp_elecs(temp_eleci)));
    does_exist = myIsField(CoI.(cur_par).(condition),elec_name);
    if does_exist == 1
        pretemp = CoI.(cur_par).(condition).data.(elec_name);
        pretemp_mask = CoI.(cur_par).(condition).sigMask.(elec_name)
        MI1 = CoI.(cur_par).(condition).MI1.(elec_name);
        MI2 = CoI.(cur_par).(condition).MI2.(elec_name);
        MI2(:,2) = [];

        if length(pretemp_mask) == 350
            pretemp(350,:) = [];
            pretemp_mask(350,:) = [];
            pretemp(:,350) = [];
            pretemp_mask(:,350) = [];
            MI1(350) =[];
            MI2(350) =[];
        end



    tempFFiall = squeeze(cat(3,tempFFiall,pretemp));
    tempFFm = tempFFm + logical(pretemp_mask);
    temp_all_mask = squeeze(cat(3,temp_all_mask,pretemp_mask));
    tempFFmr = tempFFmr + logical(pretemp_mask>0);
    tempFFms = tempFFms + logical(pretemp_mask<0);
    tempFFmi1 = squeeze(cat(2,tempFFmi1, MI1));
    tempFFmi2 = squeeze(cat(2,tempFFmi2, MI2));

    if temp_nb == 1
        tempFFmi1 = horzcat(tempFFmi1,tempFFmi1);
    end

    else
        continue
    end

end

tempFFi = squeeze(mean(tempFFiall,3));
tempFFir = tempFFi;
tempFFis = tempFFi;
tempFFir(tempFFi <= 0 ) = 0;
tempFFis(tempFFi >= 0 ) = 0;

%Frontal within
for front_eleci = 1:length(front_elecs)
    elec_name = char(strcat(front_elecs(front_eleci),'_',front_elecs(front_eleci)));
    
    does_exist = myIsField(CoI.(cur_par).(condition),elec_name);
    if does_exist == 1

        pretemp = CoI.(cur_par).(condition).data.(elec_name);
        pretemp_mask = CoI.(cur_par).(condition).sigMask.(elec_name)
        MI1 = CoI.(cur_par).(condition).MI1.(elec_name);
        MI2 = CoI.(cur_par).(condition).MI2.(elec_name)
        MI2(:,2) = [];

        if length(pretemp_mask) == 350
            pretemp(350,:) = [];
            pretemp_mask(350,:) = [];
            pretemp(:,350) = [];
            pretemp_mask(:,350) = [];
            MI1(350) =[];
            MI2(350) =[];
        end


    frontFFiall = squeeze(cat(3,frontFFiall,pretemp));
    front_all_mask = squeeze(cat(3,front_all_mask,pretemp_mask));
    frontFFm = frontFFm + logical(pretemp_mask);
    frontFFmr = frontFFmr + logical(pretemp_mask>0);
    frontFFms = frontFFms + logical(pretemp_mask<0);
    frontFFmi1 = squeeze(cat(2,frontFFmi1, MI1));
    frontFFmi2 = squeeze(cat(2,frontFFmi2, MI2));

    if front_nb == 1
        frontFFmi1 = horzcat(frontFFmi1,frontFFmi1);
    end
    else 
        continue
    end
   
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
        does_exist = myIsField(CoI.(cur_par).(condition),elec_name);
    if does_exist == 1

        pretemp = CoI.(cur_par).(condition).data.(elec_name);
        pretemp_mask = CoI.(cur_par).(condition).sigMask.(elec_name)
        MI1 = CoI.(cur_par).(condition).MI1.(elec_name);
        MI2 = CoI.(cur_par).(condition).MI2.(elec_name)
        MI2(:,2) = [];

        if length(pretemp_mask) == 350
            pretemp(350,:) = [];
            pretemp_mask(350,:) = [];
            pretemp(:,350) = [];
            pretemp_mask(:,350) = [];
            MI1(350) =[];
            MI2(350) =[];
        end


        tempfrontFFiall = squeeze(cat(3,tempfrontFFiall,pretemp));
        tempfrontFFm = tempfrontFFm + logical(pretemp_mask);
        tempfrontFFmr = tempfrontFFmr + logical(pretemp_mask>0);
        tempfrontFFms = tempfrontFFms + logical(pretemp_mask<0);
        tempfrontFFmi1 = squeeze(cat(2,tempfrontFFmi1, MI1));
        tempfrontFFmi2 = squeeze(cat(2,tempfrontFFmi2, MI2));
        nb = nb + 1;
        else 
        continue
    end
   
    end
end

tempfrontFFi = squeeze(mean(tempfrontFFiall,3));
tempfrontFFir = tempfrontFFi;
tempfrontFFis = tempfrontFFi;
tempfrontFFir(tempfrontFFi <= 0 ) = 0;
tempfrontFFis(tempfrontFFi >= 0 ) = 0;
tempfront_nb = nb;
final_tempfront = final_tempfront + tempfront_nb

%Temporo-temporal
nb = 0;
elec_num = 1;
double_cnt = 1;
pckd_cnt = 1;
all_elec_names ={};
double_elecs = {};
picked_elecs = {};

for temp_eleci = 1:length(temp_elecs)
    for temp_eleci_2 = 1:length(temp_elecs)

        %if elec1 == elec2, skip (as we dont want within elecs)
        E_one = char(temp_elecs(temp_eleci));
        E_two = char(temp_elecs(temp_eleci_2));

        if strcmp(E_one,E_two) == 1
            continue
        end

        elec_name = char(strcat(temp_elecs(temp_eleci),'_',temp_elecs(temp_eleci_2)));
        elec_name_inverse = char(strcat(temp_elecs(temp_eleci_2),'_',temp_elecs(temp_eleci)));
        all_elec_names{elec_num} = elec_name;
        if any(strcmp(all_elec_names,elec_name_inverse))
            double_elecs{double_cnt} = elec_name;
            double_cnt = double_cnt + 1;
            elec_num = elec_num + 1;
            continue

        end
        picked_elecs{pckd_cnt} = elec_name;
        does_exist = myIsField(CoI.(cur_par).(condition),elec_name);
        if does_exist == 1

        pretemp = CoI.(cur_par).(condition).data.(elec_name);
        pretemp_mask = CoI.(cur_par).(condition).sigMask.(elec_name)
        MI1 = CoI.(cur_par).(condition).MI1.(elec_name);
        MI2 = CoI.(cur_par).(condition).MI2.(elec_name)
        MI2(:,2) = [];

        if length(pretemp_mask) == 350
            pretemp(350,:) = [];
            pretemp_mask(350,:) = [];
            pretemp(:,350) = [];
            pretemp_mask(:,350) = [];
            MI1(350) =[];
            MI2(350) =[];
        end



        tempintFFiall = squeeze(cat(3,tempintFFiall,pretemp));
        tempintFFm = tempintFFm + logical(pretemp_mask);
        tempintFFmr = tempintFFmr + logical(pretemp_mask>0);
        tempintFFms = tempintFFms + logical(pretemp_mask<0);
        tempintFFmi1 = squeeze(cat(2,tempintFFmi1, MI1));
        tempintFFmi2 = squeeze(cat(2,tempintFFmi2, MI2));
        nb = nb + 1;
        elec_num = elec_num + 1;
        pckd_cnt = pckd_cnt + 1;
         else 
        continue
    end
    end
end

tempintFFi = squeeze(mean(tempintFFiall,3));
tempintFFir = tempintFFi;
tempintFFis = tempintFFi;
tempintFFir(tempintFFi <= 0 ) = 0;
tempintFFis(tempintFFi >= 0 ) = 0;
tempint_nb = nb
final_tempint = final_tempint + tempint_nb

%Fronto-frontal
nb = 0;
elec_num = 1;
double_cnt = 1;
pckd_cnt = 1;
all_elec_names ={};
double_elecs = {};
picked_elecs = {};
for front_eleci = 1:length(front_elecs)
    for front_eleci_2 = 1:length(front_elecs)

        %if elec1 == elec2, skip (as we dont want within elecs)
        E_one = char(front_elecs(front_eleci));
        E_two = char(front_elecs(front_eleci_2));

        if strcmp(E_one,E_two) == 1
            continue
        end


        elec_name = char(strcat(front_elecs(front_eleci),'_',front_elecs(front_eleci_2)));
        elec_name_inverse = char(strcat(front_elecs(front_eleci_2),'_',front_elecs(front_eleci)));
        all_elec_names{elec_num} = elec_name;
        if any(strcmp(all_elec_names,elec_name_inverse))
            double_elecs{double_cnt} = elec_name;
            double_cnt = double_cnt + 1;
            elec_num = elec_num + 1;
            continue

        end

        picked_elecs{pckd_cnt} = elec_name;
        does_exist = myIsField(CoI.(cur_par).(condition),elec_name);
        if does_exist == 1

        pretemp = CoI.(cur_par).(condition).data.(elec_name);
        pretemp_mask = CoI.(cur_par).(condition).sigMask.(elec_name)
        MI1 = CoI.(cur_par).(condition).MI1.(elec_name);
        MI2 = CoI.(cur_par).(condition).MI2.(elec_name)
        MI2(:,2) = [];

        if length(pretemp_mask) == 350
            pretemp(350,:) = [];
            pretemp_mask(350,:) = [];
            pretemp(:,350) = [];
            pretemp_mask(:,350) = [];
            MI1(350) =[];
            MI2(350) =[];
        end


        frontintFFiall = squeeze(cat(3,frontintFFiall,pretemp));
        frontintFFm = frontintFFm + logical(pretemp_mask);
        frontintFFmr = frontintFFmr + logical(pretemp_mask>0);
        frontintFFms = frontintFFms + logical(pretemp_mask<0);
        frontintFFmi1 = squeeze(cat(2,frontintFFmi1, MI1));
        frontintFFmi2 = squeeze(cat(2,frontintFFmi2, MI2));
        nb = nb + 1;
        elec_num = elec_num + 1;
        pckd_cnt = pckd_cnt + 1;
         else 
        continue
    end
    end
end
all_elec_names;
frontintFFi = squeeze(mean(frontintFFiall,3));
frontintFFir = frontintFFi;
frontintFFis = frontintFFi;
frontintFFir(frontintFFi <= 0 ) = 0;
frontintFFis(frontintFFi >= 0 ) = 0;
frontint_nb = nb
final_frontint = final_frontint + frontint_nb;

end
end





