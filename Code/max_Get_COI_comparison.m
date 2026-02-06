%% GET CO-INFORMATION
% Estimates the Co-information for the specified electrode pair via
% Gaussian copula estimation (Ince et al., 2016)


function [MI_stat]  = max_Get_COI_comparison(electrodes, basefold, datatype,...
    subject, condition, deviant_group_number_1, ...
    standard_group_number_1, deviant_group_number_2, standard_group_number_2,...
     corrected, srate, baseline, kperm, threshold)

overVar_file = strcat(basefold, datatype, '/', subject, '_', condition, '.mat');
[dvt_1, std_1] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_1, standard_group_number_1,...
    corrected, srate); %Load deviant and standard trials for activity 1
[dvt_2, std_2] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_2, standard_group_number_2,...
    corrected, srate); %Load deviant and standard trials for activity 2

patname = char(subject);

%% CHECK FOR EQUAL CHANNELS & TRIALS

if (std_1.trials < std_2.trials) %If there are fewer activity 1 trials, reduce activity 2 trials
    trialnum = std_1.trials;
    std_2.data = std_2.data(:,:,1:trialnum);
    std_2.trials = trialnum;
    dvt_2.data = dvt_2.data(:,:,1:trialnum);
    dvt_2.trials = trialnum;
else %Otherwise do the opposite
    trialnum = std_2.trials;
    std_1.data = std_1.data(:,:,1:trialnum);
    std_1.trials = trialnum;
    dvt_1.data = dvt_1.data(:,:,1:trialnum);
    dvt_1.trials = trialnum;
end

%% GET DATA & BASELINE NORMALISE


%First normalise for activity 1
data_dvt_1 = permute(dvt_1.data,[1 3 2]); mean_data_dvt_1 = squeeze(mean(data_dvt_1(:,:,baseline),3));
dvt_1.data = data_dvt_1 - repmat(mean_data_dvt_1,[1 1 size(dvt_1.data,2)]);
bb1_dev_bl_1 = permute(dvt_1.data,[1 3 2]);

data_std_1 = permute(std_1.data,[1 3 2]); mean_data_std_1 = squeeze(mean(data_std_1(:,:,baseline),3));
std_1.data = data_std_1 - repmat(mean_data_std_1,[1 1 size(std_1.data,2)]);
bb1_std_bl_1 = permute(std_1.data,[1 3 2]);

dvt_1_data_minus_mean = permute(bb1_dev_bl_1,[3 1 2]);
std_1_data_minus_mean = permute(bb1_std_bl_1,[3 1 2]);

% Then normalise for activity 2
data_dvt_2 = permute(dvt_2.data,[1 3 2]); mean_data_dvt_2 = squeeze(mean(data_dvt_2(:,:,baseline),3));
dvt_2.data = data_dvt_2 - repmat(mean_data_dvt_2,[1 1 size(dvt_2.data,2)]);
bb1_dev_bl_2 = permute(dvt_2.data,[1 3 2]);

data_std_2 = permute(std_2.data,[1 3 2]); mean_data_std_2 = squeeze(mean(data_std_2(:,:,baseline),3));
std_2.data = data_std_2 - repmat(mean_data_std_2,[1 1 size(std_2.data,2)]);
bb1_std_bl_2 = permute(std_2.data,[1 3 2]);

dvt_2_data_minus_mean = permute(bb1_dev_bl_2,[3 1 2]);
std_2_data_minus_mean = permute(bb1_std_bl_2,[3 1 2]);


%% GET THE RIGHT CHANNELS
% looks up the data for thes specified channels
if isempty(std_1.chanlocs) == 1
    all_elecs ={};
    for e_index = 1:std_1.nbchan
        name = strcat('E',num2str(e_index));
        all_elecs(e_index) = {name};
    end
    All_elecs =  all_elecs';
else
    All_elecs = std_1.chanlocs.labels;
end

channel1 = extractBefore(electrodes,'_permuted_');
channel2 = extractAfter (electrodes,'_permuted_');
chnum1 = extractBefore(electrodes,'_permuted_');
chnum2 = extractAfter (electrodes,'_permuted_');
channel1 = strcmp(channel1, All_elecs');
channel2 = strcmp(channel2, All_elecs');
comb_chan = char(strcat(chnum1,'_',chnum2));
comb_chan = replace(comb_chan,'-','_');
%Get correct channel data for activity tag 1
dvt_1_data_ERP_sel_chan1 = squeeze(dvt_1_data_minus_mean(:, channel1,:));
std_1_data_ERP_sel_chan1 = squeeze(std_1_data_minus_mean(:,channel1,:));
dvt_1_data_ERP_sel_chan2 = squeeze(dvt_1_data_minus_mean(:, channel2, :));
std_1_data_ERP_sel_chan2 = squeeze(std_1_data_minus_mean(:,channel2,:));
%Get correct channel data for activity tag 2
dvt_2_data_ERP_sel_chan1 = squeeze(dvt_2_data_minus_mean(:, channel1,:));
std_2_data_ERP_sel_chan1 = squeeze(std_2_data_minus_mean(:,channel1,:));
dvt_2_data_ERP_sel_chan2 = squeeze(dvt_2_data_minus_mean(:, channel2, :));
std_2_data_ERP_sel_chan2 = squeeze(std_2_data_minus_mean(:,channel2,:));

%% COMPUTE CO-INFORMATION WITH GCMI
% Computes co-information for the specified channels via gaussian copula
% estimation (Ince et al., 2016)

%Computes the CoI of channel X against channel X (within)
if strcmp(chnum1,chnum2) == 1
    % [MI1, MI2, coI, jointMI, sigMask] = cnm_coI_stimtimetime_within_channel([dvt_data_ERP_sel_chan1; std_data_ERP_sel_chan1],[dvt_data_ERP_sel_chan1; std_data_ERP_sel_chan1], [zeros(1, size(dvt_data_minus_mean,1)), ones(1, size(dvt_data_minus_mean,1))]', kperm);
    [MI1_diff, MI2_diff, diffCoI, diffJointMI, sigMask] = cnm_coI_stimtimetime_perm_diff_within_channel(...
     [dvt_1_data_ERP_sel_chan1; std_1_data_ERP_sel_chan1],...
     [dvt_1_data_ERP_sel_chan2; std_1_data_ERP_sel_chan2],...
     [dvt_2_data_ERP_sel_chan1; std_2_data_ERP_sel_chan1],...
     [dvt_2_data_ERP_sel_chan2; std_2_data_ERP_sel_chan2],...
     [zeros(1, size(dvt_1_data_minus_mean,1)), ones(1, size(dvt_1_data_minus_mean,1))]', kperm, threshold);

    %Computes the CoI of channel X against channel Y (between)
else
    % [MI1, MI2, coI, jointMI, sigMask] = cnm_coI_stimtimetime([dvt_data_ERP_sel_chan1; std_data_ERP_sel_chan1],[dvt_data_ERP_sel_chan2; std_data_ERP_sel_chan2], [zeros(1, size(dvt_data_minus_mean,1)), ones(1, size(dvt_data_minus_mean,1))]', kperm);
    [MI1_diff, MI2_diff, diffCoI, diffJointMI, sigMask] = cnm_coI_stimtimetime_perm_diff(...
     [dvt_1_data_ERP_sel_chan1; std_1_data_ERP_sel_chan1],...
     [dvt_1_data_ERP_sel_chan2; std_1_data_ERP_sel_chan2],...
     [dvt_2_data_ERP_sel_chan1; std_2_data_ERP_sel_chan1],...
     [dvt_2_data_ERP_sel_chan2; std_2_data_ERP_sel_chan2],...
     [zeros(1, size(dvt_1_data_minus_mean,1)), ones(1, size(dvt_1_data_minus_mean,1))]', kperm, threshold);

end
MI_stat.CoI.(patname).(char(condition)).diffCoI.(comb_chan) = diffCoI;
MI_stat.CoI.(patname).electrode.(char(condition)).(comb_chan) = electrodes;
MI_stat.CoI.(patname).(char(condition)).sigMask.(comb_chan)= sigMask;
MI_stat.CoI.(patname).(char(condition)).MI1_diff.(comb_chan)= MI1_diff;
MI_stat.CoI.(patname).(char(condition)).MI2_diff.(comb_chan)= MI2_diff;
MI_stat.CoI.(patname).(char(condition)).jointJointMI.(comb_chan)= diffJointMI;
end