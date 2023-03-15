%Gets ERPs for specified condition and participant

[basefold, datatype, participantnum, ~, condition, participants, ~, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(0);
[std, dvt] = impiEEG(participantnum, basefold, datatype, condition, srate, low_cutoff, high_cutoff, filt_order,re_epoch, dev_epochs, std_epochs, epoch_length);


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

        
%% PLOT ERPS

        dvt_data_minus_mean = permute(dvt_data_minus_mean,[2 3 1]);
        std_data_minus_mean = permute(std_data_minus_mean,[2 3 1]);
        H = figure(1);
        for ch = 1 : dvt.nbchan
           nexttile(ch)
           stdshade_acj(squeeze(std_data_minus_mean(ch,:,:))',0.2,'k',dvt.times); hold on
           stdshade_acj(squeeze(dvt_data_minus_mean(ch,:,:))',0.2,'g',dvt.times);
            set(gca,'xlim',[700 1150],'ylim',[-100 100],'xtick',[800 950 1100], 'xticklabel', [0, 150, 300],'ytick',[-100 100], 'yticklabel', [-100 100]);
           hold off
        end
%________________________________________________________________________________________________________________________________________________________
