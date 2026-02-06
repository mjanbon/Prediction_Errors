%% ESTIMATES MI AND GETS ELECTRODES OF INTEREST
% Estimates MI via gaussian copula estimation (Ince et al., 2016), performs
% parametric permutation testing and saves electrodes with significant MI
% into a mat-file.


function main_MI_ERP_comparison(task_id)
USING_HPC = 2; % 0 for local, 1 for Cambridge HPC, 2 for QMUL HPC
%Get the parameters for the participant and datatype you want
% [basefold, datatype, subject, all_con, condition, participants, ~ , re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(0);
if USING_HPC == 1
    addpath(genpath('/home/mj649/rds/hpc-work/CNM')); % Add matlab paths to code folders and subfolders
    addpath(genpath('/home/mj649/rds/hpc-work/GCMI_master'));
    addpath(genpath('/home/mj649/rds/hpc-work/Prediction_Errors'));
end

if USING_HPC == 2 % Add matlab paths to code folders and subfolders on QMUL Apocrita
    addpath(genpath('/data/SBBS-PIDProject/Maxime/CNM')); % Add matlab paths to code folders and subfolders
    addpath(genpath('/data/SBBS-PIDProject/Maxime/GCMI_master'));
    addpath(genpath('/data/SBBS-PIDProject/Maxime/Prediction_Errors'));
end
THRESHOLD = 0.05
[basefold, datatype,all_con, condition,subject,participants, EoI,...
    srate, activity_tag_1, deviant_group_number_1, standard_group_number_1,...
    activity_tag_2, deviant_group_number_2, standard_group_number_2, corrected,...
    stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_comparison_param(USING_HPC, 0)



%%


%Go through all specified participants and conditions
for i = 1: length(participants)
    for con = 1:length(all_con)
        %% IMPORT DATA
        participants(i)
        all_con(con)
        % Import data and check for same amount of trials & channels
        overVar_file = strcat(basefold, datatype, '/', participants(i), '_', all_con(con), '.mat')
        if exist(char(overVar_file)) == 0
            continue
        end
        [dvt_1, std_1] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_1, standard_group_number_1,...
            corrected, srate); %Load deviant and standard trials for activity 1
        [dvt_2, std_2] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_2, standard_group_number_2,...
            corrected, srate); %Load deviant and standard trials for activity 2

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


        %% BASELINE NORMALISATION
        %First normalise for activity 1
        data_dvt_1 = permute(dvt_1.data,[1 3 2]); mean_data_dvt_1 = squeeze(mean(data_dvt_1(:,:,baseline),3));
        dvt_1.data = data_dvt_1 - repmat(mean_data_dvt_1,[1 1 size(dvt_1.data,2)]);
        bb1_dev_bl_1 = permute(dvt_1.data,[1 3 2]);

        data_std_1 = permute(std_1.data,[1 3 2]); mean_data_std_1 = squeeze(mean(data_std_1(:,:,baseline),3));
        std_1.data = data_std_1 - repmat(mean_data_std_1,[1 1 size(std_1.data,2)]);
        bb1_std_bl_1 = permute(std_1.data,[1 3 2]);

        bb_dev_1 = permute(bb1_dev_bl_1,[3 1 2]);
        bb_std_1 = permute(bb1_std_bl_1,[3 1 2]);

        % Then normalise for activity 2
        data_dvt_2 = permute(dvt_2.data,[1 3 2]); mean_data_dvt_2 = squeeze(mean(data_dvt_2(:,:,baseline),3));
        dvt_2.data = data_dvt_2 - repmat(mean_data_dvt_2,[1 1 size(dvt_2.data,2)]);
        bb1_dev_bl_2 = permute(dvt_2.data,[1 3 2]);

        data_std_2 = permute(std_2.data,[1 3 2]); mean_data_std_2 = squeeze(mean(data_std_2(:,:,baseline),3));
        std_2.data = data_std_2 - repmat(mean_data_std_2,[1 1 size(std_2.data,2)]);
        bb1_std_bl_2 = permute(std_2.data,[1 3 2]);

        bb_dev_2 = permute(bb1_dev_bl_2,[3 1 2]);
        bb_std_2 = permute(bb1_std_bl_2,[3 1 2]);

        %% CALCULATE MUTUAL INFORMATION & PERFORM PERMUTATION TESTING USING GCMI

        elec_of_I = {};
        h = figure(con);
        for ch = 1 : dvt_1.nbchan %Still go through every channel
            chnum = strcat('ch_',string(ch));
            dat_1.class1.BB = squeeze(bb_dev_1(:, ch, :)); %Deviant trials for activity 1 channel ch
            dat_1.class2.BB = squeeze(bb_std_1(:, ch, :)); %Standard trials for activity 1 channel ch

            dat_2.class1.BB = squeeze(bb_dev_2(:, ch, :)); %Deviant trials for activity 2 channel ch
            dat_2.class2.BB = squeeze(bb_std_2(:, ch, :)); %Standard trials for activity 2 channel ch

            %Calculate the difference in MI between the specified channels
            [MI_diff, MI_1, MI_2, sigMask] = cnm_MI_stimtime_perm_diff([dat_1.class1.BB; dat_1.class2.BB],[dat_2.class1.BB; dat_2.class2.BB],...
                [zeros(1, size(bb_dev_1,1)), ones(1, size(bb_dev_1,1))]', kperm, THRESHOLD);

            %save
            chan_name = strcat('E',num2str(ch));
            MI_stat.MI_diff.(chan_name) = MI_diff;
            MI_stat.MI_1.(chan_name) = MI_1;
            MI_stat.MI_2.(chan_name) = MI_2;
            MI_stat.electrode = chan_name;
            MI_stat.sigMask.(chan_name)= sigMask;

            %% PLOT
            %get correct timing
            timing = dvt_1.times;
            % timing(start_cut_off) = [];
            % timing(end_cut_off:end) = [];
            % timing = timing';

            %plot MI
            % nexttile(ch)
            % plot(timing, MI_stat.(participants{i}).(all_con{con}).MI.(chan_name));
            % title(chan_name)
            % hold on;
            % xlabel('Time [ms]');
            % yLimits = [-0.001 0.2];
            % pos_sigbar = yLimits(2) - (0.07 * range(yLimits));
            % xline(stim_onset, '--');


            % % Draw stats
            % stat= logical(MI_stat.sigMask.(chan_name));
            % ylim(yLimits);
            % ylabel('Mutual information (bits)');

            EoI_list = true;
            % for tIdx = 1:length(timing)-1
            %     tIdx2 = timing(tIdx);
            %     tIdx2_1 = timing(tIdx+1);
            %     if stat(tIdx) > 0
            %         plot([tIdx2, tIdx2_1], [pos_sigbar, pos_sigbar], 'LineWidth', 3, 'Color', 'm');
            %         EoI_list = true;
            %     end
            % end

            %List electrodes of interest
            if EoI_list == true
                % chan = (MI_stat.(participants{i}).(all_con{con}).electrode); %Not electrode(i) anymore
                % chan = replace(chan,'-','_');
                Electorodes.(chan_name) = chan_name;
            end
            hold off
        end
        %
        %% PLOT ERPS

        % H = figure(length(all_con)+1);
        % for ch = 1 : dvt.nbchan
        %    nexttile(ch)
        %    chan_name = strcat('E',num2str(ch));
        %    plot(dvt.times, squeeze(mean(dvt.data(ch,:,:),2)), 'DisplayName', 'D'); hold on
        %    %stdshade_acj(squeeze(dvt_data_minus_mean(ch,:,:))',0.2,'g',dvt.times);
        %    plot(std.times, squeeze(mean(std.data(ch,:,:),2)), 'DisplayName','S'); hold on
        %    title(chan_name)
        %    xlabel('Time [ms]');
        %    ylabel('LFP');
        %    xline(stim_onset, '--');
        %    legend('S','D')
        %    hold off
        % end
        %% SAVE ALL
        % basename = strcat (participants(i), 'MIs.mat');
        % EoI_names = strcat (participants(i), 'EoI.mat');
        % cd (basefold)
        % c_tag='NC';
        % if corrected==1
        %     c_tag='C';
        % end
        % cd('MI_Figures') % Save MI figures for the fly
        % set(h, 'Renderer', 'painters');
        % print(h, strcat(char(participants(i)),'_', char(activity_tag), '_',char(all_con(con)),...
        %     c_tag, num2str(standard_group_number),...
        %     num2str(deviant_group_number),'_MI_figures_local'), '-dsvg');
        % saveas(H,strcat(char(participants(i)),'_', char(activity_tag), '_', char(all_con(con)),...
        %     c_tag, num2str(standard_group_number),...
        %     num2str(deviant_group_number),'_ERPs_local'),'svg');
        % cd ../

        MI_name = char(strcat(participants(i),(strcat(char(activity_tag_1),'_',char(activity_tag_2))),char(all_con(con)),'_MI_data.mat'));
        save (fullfile(basefold, 'MI_Data', MI_name),'MI_stat','-mat')

        
        % cd('DataEoI') % Save EoI data for the fly
        filename = fullfile(basefold, 'DataEoI', ['EoI_data_' datatype '.mat']);

        if exist(filename, 'file')
            load(filename, 'EoI');
        else
            EoI = struct();
        end

        EoI.(char(participants(i))).(strcat(char(activity_tag_1),'_',char(activity_tag_2))).(char(all_con(con))) = Electorodes;
        EoI.(char(participants(i))).(strcat(char(activity_tag_1),'_',char(activity_tag_2))).(char(all_con(con))) = fieldnames(EoI.(char(participants(i))).(strcat(char(activity_tag_1),'_',char(activity_tag_2))).(char(all_con(con))));
        save(filename, 'EoI', '-mat');

        clear Electorodes
        clear MI_stat
        close all

    end
end
end




