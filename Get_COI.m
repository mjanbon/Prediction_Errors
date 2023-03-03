%% GET CO-INFORMATION
% Estimates the Co-information for the specified electrode pair via
% Gaussian copula estimation (Ince et al., 2016)


function [MI_stat]  = Get_COI(electrodes,basefold, datatype, participant, condition, participants, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm)
addpath('/home/jma201/coi-ieeg')
addpath('/home/jma201/iEEG')

[std, dvt] = impiEEG(participant, basefold, datatype, condition, srate, low_cutoff, high_cutoff, filt_order, re_epoch, dev_epochs, std_epochs, epoch_length);

patname = char(participants(participant));

%% CHECK FOR EQUAL CHANNELS & TRIALS

if isempty(dvt) == 1
    fprintf('dvt is empty')
    return

elseif (dvt.nbchan < std.nbchan)
    channum = dvt.nbchan;
    std = pop_select(std, 'channel', 1:channum);
else
    channum = std.nbchan;
    dvt = pop_select(dvt, 'channel', 1:channum);
end

if isempty(dvt) == 1
    fprintf('dvt is empty')
    return

elseif (dvt.trials < std.trials)
    trialnum = dvt.trials;
    std = pop_select(std, 'trial', 1:trialnum);
else
    trialnum = std.trials;
    dvt = pop_select(dvt, 'trial', 1:trialnum);
end

if isempty(std) == 1
    fprintf('dvt is empty')
    return

elseif (std.nbchan < dvt.nbchan)
    channum = std.nbchan;
    dvt = pop_select(dvt, 'channel', 1:channum);
else
    channum = dvt.nbchan;
    std = pop_select(std, 'channel', 1:channum);
end

if isempty(dvt) == 1
    fprintf('dvt is empty')
    return

elseif (std.trials < dvt.trials)
    trialnum = std.trials;
    dvt = pop_select(dvt, 'trial', 1:trialnum);
else
    trialnum = std.trials;
    std = pop_select(std, 'trial', 1:trialnum);
end

%% GET DATA & BASELINE NORMALISE
%Get time series data (electrode x time x trial) into a matrix
std_data = std.data;
dvt_data = dvt.data;

%Baseline normalisation for std_data
std_reorder = permute(std_data,[1, 3, 2]);
mean_std_trials = squeeze(mean(std_reorder(:,:,baseline),3));
mean_std_trials_repmat_dimension = repmat(mean_std_trials,[1 1 size(std_data,2)]);
std_data_minus_mean = std_reorder - mean_std_trials_repmat_dimension;
std_data_minus_mean = permute(std_data_minus_mean, [ 1 3 2]);

%Baseline normalisation for dvt_data
dvt_reorder = permute(dvt_data,[1, 3, 2]);
mean_dvt_trials = squeeze(mean(dvt_reorder(:,:,baseline),3));
mean_dvt_trials_repmat_dimension = repmat(mean_dvt_trials,[1 1 size(dvt_data,2)]);
dvt_data_minus_mean = dvt_reorder - mean_dvt_trials_repmat_dimension;
dvt_data_minus_mean = permute(dvt_data_minus_mean, [ 1 3 2]);

%Permute data to trials x electrodes x timepoints
dvt_data_minus_mean = permute(dvt_data_minus_mean,[3 1 2]);
std_data_minus_mean = permute(std_data_minus_mean,[3 1 2]);

%Get rid of timepoints we are not interested in for efficiency
dvt_data_minus_mean(:,:,start_cut_off) = [];
std_data_minus_mean(:,:,start_cut_off) = [];
dvt.times (start_cut_off) = [];

dvt_data_minus_mean(:,:,end_cut_off:end) = [];
std_data_minus_mean(:,:,end_cut_off:end) = [];
dvt.times (end_cut_off:end) = [];


%% GET THE RIGHT CHANNELS
% looks up the data for thes specified channels
if isempty(std.chanlocs) == 1
    all_elecs ={};
    for e_index = 1:std.nbchan
        name = strcat('E',num2str(e_index));
        all_elecs(e_index) = {name};
    end
    All_elecs =  all_elecs';
else
    All_elecs = std.chanlocs.labels;
end

channel1 = extractBefore(electrodes,'_permuted_');
channel2 = extractAfter (electrodes,'_permuted_');
chnum1 = extractBefore(electrodes,'_permuted_');
chnum2 = extractAfter (electrodes,'_permuted_');
channel1 = strcmp(channel1, All_elecs');
channel2 = strcmp(channel2, All_elecs');
comb_chan = char(strcat(chnum1,'_',chnum2));
comb_chan = replace(comb_chan,'-','_');
dvt_data_ERP_sel_chan1 = squeeze(dvt_data_minus_mean(:, channel1,:));
std_data_ERP_sel_chan1 = squeeze(std_data_minus_mean(:,channel1,:));
dvt_data_ERP_sel_chan2 = squeeze(dvt_data_minus_mean(:, channel2, :));
std_data_ERP_sel_chan2 = squeeze(std_data_minus_mean(:,channel2,:));

%% COMPUTE CO-INFORMATION WITH GCMI
% Computes co-information for the specified channels via gaussian copula
% estimation (Ince et al., 2016)

if strcmp(chnum1,chnum2) == 1
    [MI1, MI2, coI, jointMI, sigMask] = cnm_coI_stimtimetime_within_channel([dvt_data_ERP_sel_chan1; std_data_ERP_sel_chan1],[dvt_data_ERP_sel_chan1; std_data_ERP_sel_chan1], [zeros(1, size(dvt_data_minus_mean,1)), ones(1, size(dvt_data_minus_mean,1))]', kperm);

    MI_stat.CoI.(patname).(char(condition)).data.(comb_chan) = coI;
    MI_stat.CoI.(patname).electrode.(char(condition)).(comb_chan) = electrodes;
    MI_stat.CoI.(patname).(char(condition)).sigMask.(comb_chan)= sigMask;
    MI_stat.CoI.(patname).(char(condition)).MI1.(comb_chan)= MI1;
    MI_stat.CoI.(patname).(char(condition)).MI2.(comb_chan)= MI2;
    MI_stat.CoI.(patname).(char(condition)).joint.(comb_chan)= jointMI;

    %Computes the CoI of channel X against channel Y (between)
else
    [MI1, MI2, coI, jointMI, sigMask] = cnm_coI_stimtimetime([dvt_data_ERP_sel_chan1; std_data_ERP_sel_chan1],[dvt_data_ERP_sel_chan2; std_data_ERP_sel_chan2], [zeros(1, size(dvt_data_minus_mean,1)), ones(1, size(dvt_data_minus_mean,1))]', kperm);

    MI_stat.CoI.(patname).(char(condition)).data.(comb_chan) = coI;
    MI_stat.CoI.(patname).electrode.(char(condition)).(comb_chan) = electrodes;
    MI_stat.CoI.(patname).(char(condition)).sigMask.(comb_chan)= sigMask;
    MI_stat.CoI.(patname).(char(condition)).MI1.(comb_chan)= MI1;
    MI_stat.CoI.(patname).(char(condition)).MI2.(comb_chan)= MI2;
    MI_stat.CoI.(patname).(char(condition)).joint.(comb_chan)= jointMI;

end

end



