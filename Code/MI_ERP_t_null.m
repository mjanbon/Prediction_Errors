%% ESTIMATES MI AND GETS ELECTRODES OF INTEREST
% Estimates MI via gaussian copula estimation (Ince et al., 2016), performs
% parametric permutation testing and saves electrodes with significant MI
% into a mat-file.

clear all
eeglab

%Get the parameters for the participant and datatype you wants
[basefold, datatype, subject, all_con, condition, participants, ~ , re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(0);

Sim_MI_all_par = zeros(length(participants),225);
%Go through all specified participants and conditions
for i = 1: length(participants)
    for con = 1
        %% IMPORT DATA
        participants(i)
        all_con(con)
        % Import data and check for same amount of trials & channels
        [dvt, std] = impiEEG(i, basefold, datatype, all_con(con), srate, low_cutoff, high_cutoff, filt_order,re_epoch, dev_epochs, std_epochs, epoch_length);

        if (dvt.nbchan < std.nbchan)
            channum = dvt.nbchan;
            std = pop_select(std, 'channel', 1:channum);
        else
            channum = std.nbchan;
            dvt = pop_select(dvt, 'channel', 1:channum);
        end

        if (dvt.trials < std.trials)
            trialnum = dvt.trials;
            std = pop_select(std, 'trial', 1:trialnum);
        else
            trialnum = std.trials;
            dvt = pop_select(dvt, 'trial', 1:trialnum);
        end

        if (std.nbchan < dvt.nbchan)
            channum = std.nbchan;
            dvt = pop_select(dvt, 'channel', 1:channum);
        else
            channum = dvt.nbchan;
            std = pop_select(std, 'channel', 1:channum);
        end

        if (std.trials < dvt.trials)
            trialnum = std.trials;
            dvt = pop_select(dvt, 'trial', 1:trialnum)
        else
            trialnum = std.trials;
            std = pop_select(std, 'trial', 1:trialnum)
        end

        %% BASELINE NORMALISATION
        data_dvt = permute(dvt.data,[1 3 2]); mean_data_dvt = squeeze(mean(data_dvt(:,:,baseline),3));
        dvt.data = data_dvt - repmat(mean_data_dvt,[1 1 size(dvt.data,2)]);
        dvt.data(:,:,start_cut_off) = [];
        dvt.data(:,:,end_cut_off:end) = [];
        bb1_dev_bl = permute(dvt.data,[1 3 2]);

        data_std = permute(std.data,[1 3 2]); mean_data_std = squeeze(mean(data_std(:,:,baseline),3));
        std.data = data_std - repmat(mean_data_std,[1 1 size(std.data,2)]);
        std.data(:,:,start_cut_off) = [];
        std.data(:,:,end_cut_off:end) = [];
        bb1_std_bl = permute(std.data,[1 3 2]);

        bb_dev = permute(bb1_dev_bl,[3 1 2]);
        bb_std = permute(bb1_std_bl,[3 1 2]);

    %% CALCULATE MUTUAL INFORMATION & PERFORM PERMUTATION TESTING USING GCMI
    
    elec_of_I = {};
    PermMI = zeros(dvt.nbchan, length(bb_std));
    all_mean_MI = zeros (10,length(bb_std));
    for perm =  1 : 10
    for ch = 1 : dvt.nbchan
        chnum = strcat('ch_',string(ch));
        dat.class1.BB = squeeze(bb_dev(:, ch, :));
        dat.class2.BB = squeeze(bb_std(:, ch, :));

        %Calculate MI of the specified channels (1 permutation)
        [MI_perm] = cnm_MI_stimtime_t_null([dat.class1.BB; dat.class2.BB],[zeros(1, size(bb_dev,1)), ones(1, size(bb_dev,1))]', 1);
        PermMI(ch,:) = MI_perm;
        PermMI_mean = mean(PermMI);
     
    end
    all_mean_MI(perm,:) = PermMI_mean;  
    end
    end
    Sim_MI_all_par(i,:) = all_mean_MI(randi(10),:);
end
filename = strcat('EoI_data_',datatype,'.mat');
save (filename,'EoI','-mat')




