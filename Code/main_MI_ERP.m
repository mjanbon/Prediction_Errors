%% ESTIMATES MI AND GETS ELECTRODES OF INTEREST
% Estimates MI via gaussian copula estimation (Ince et al., 2016), performs
% parametric permutation testing and saves electrodes with significant MI
% into a mat-file.


function main_MI_ERP(task_id)
USING_HPC = 1;
%Get the parameters for the participant and datatype you want
% [basefold, datatype, subject, all_con, condition, participants, ~ , re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, baseline, start_cut_off, end_cut_off, kperm] = Get_param(0);

[basefold, datatype, all_con, condition, subject,participants, EoI,...
    srate, activity_tag, deviant_group_number, standard_group_number, corrected, stim_onset, baseline,...
    start_cut_off, end_cut_off, kperm] = Max_get_param(USING_HPC, 0);


%%


%Go through all specified participants and conditions
for i = 1: length(participants)
    for con = 1:length(all_con)
        %% IMPORT DATA
        participants(i)
        all_con(con)
        % Import data and check for same amount of trials & channels
        % [dvt, std] = impiEEG(i, basefold, datatype, all_con(con), srate, low_cutoff, high_cutoff, filt_order,re_epoch, dev_epochs, std_epochs, epoch_length);
        overVar_file = strcat(basefold, datatype, '/', participants(i), '_', all_con(con), '.mat')
        [dvt, std] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number, standard_group_number,...
            corrected, srate);
        % if (dvt.nbchan < std.nbchan)
        %     channum = dvt.nbchan;
        %     std = pop_select(std, 'channel', 1:channum);
        % else
        %     channum = std.nbchan;
        %     dvt = pop_select(dvt, 'channel', 1:channum);
        % end

        if (dvt.trials < std.trials)
            trialnum = dvt.trials;
            % std = pop_select(std, 'trial', 1:trialnum);
            std.data = std.data(:,:,1:trialnum);
        else
            trialnum = std.trials;
            % dvt = pop_select(dvt, 'trial', 1:trialnum);
            dvt.data = dvt.data(:,:,1:trialnum);
        end


        %% BASELINE NORMALISATION
        data_dvt = permute(dvt.data,[1 3 2]); mean_data_dvt = squeeze(mean(data_dvt(:,:,baseline),3));
        dvt.data = data_dvt - repmat(mean_data_dvt,[1 1 size(dvt.data,2)]);
        % dvt.data(:,:,start_cut_off) = [];
        % dvt.data(:,:,end_cut_off:end) = [];
        bb1_dev_bl = permute(dvt.data,[1 3 2]);

        data_std = permute(std.data,[1 3 2]); mean_data_std = squeeze(mean(data_std(:,:,baseline),3));
        std.data = data_std - repmat(mean_data_std,[1 1 size(std.data,2)]);
        % std.data(:,:,start_cut_off) = [];
        % std.data(:,:,end_cut_off:end) = [];
        bb1_std_bl = permute(std.data,[1 3 2]);

        bb_dev = permute(bb1_dev_bl,[3 1 2]);
        bb_std = permute(bb1_std_bl,[3 1 2]);

        %% CALCULATE MUTUAL INFORMATION & PERFORM PERMUTATION TESTING USING GCMI

        elec_of_I = {};
        h = figure(con);
        for ch = 1 : dvt.nbchan
            chnum = strcat('ch_',string(ch));
            dat.class1.BB = squeeze(bb_dev(:, ch, :));
            dat.class2.BB = squeeze(bb_std(:, ch, :));

            %Calculate MI between the specified channels
            [MI, sigMask] = cnm_MI_stimtime([dat.class1.BB; dat.class2.BB],[zeros(1, size(bb_dev,1)), ones(1, size(bb_dev,1))]', kperm);

            %save
            chan_name = strcat('E',num2str(ch));
            MI_stat.(participants{i}).(all_con{con}).MI.(chan_name) = MI;
            MI_stat.(participants{i}).(all_con{con}).electrode = chan_name;
            MI_stat.(participants{i}).(all_con{con}).sigMask.(chan_name)= sigMask;

            %% PLOT
            %get correct timing
            timing = dvt.times;
            % timing(start_cut_off) = [];
            % timing(end_cut_off:end) = [];
            % timing = timing';

            %plot MI
            nexttile(ch)
            plot(timing, MI_stat.(participants{i}).(all_con{con}).MI.(chan_name));
            title(chan_name)
            hold on;
            xlabel('Time [ms]');
            yLimits = [-0.001 0.2];
            %         xlim([0.7 1.150])
            pos_sigbar = yLimits(2) - (0.07 * range(yLimits));
            xline(stim_onset, '--');


            % Draw stats
            stat= logical(MI_stat.(participants{i}).(all_con{con}).sigMask.(chan_name));
            ylim(yLimits);
            ylabel('Mutual information (bits)');

            EoI_list = false;
            for tIdx = 1:length(timing)-1
                tIdx2 = timing(tIdx);
                tIdx2_1 = timing(tIdx+1);
                if stat(tIdx) > 0
                    plot([tIdx2, tIdx2_1], [pos_sigbar, pos_sigbar], 'LineWidth', 3, 'Color', 'm');
                    EoI_list = true;
                end
            end

            %List electrodes of interest
            if EoI_list == true
                chan = (MI_stat.(participants{i}).(all_con{con}).electrode(i));
                chan = replace(chan,'-','_');
                Electorodes.(chan_name) = chan_name;
            end
            hold off
        end
        %
        %% PLOT ERPS

        % dvt_data_minus_mean = permute(dvt_data_minus_mean,[2 3 1]);
        % std_data_minus_mean = permute(std_data_minus_mean,[2 3 1]);
        H = figure(2);
        for ch = 1 : dvt.nbchan
           nexttile(ch)
           chan_name = strcat('E',num2str(ch));
           plot(dvt.times, squeeze(mean(dvt.data(ch,:,:),2)), 'DisplayName', 'D'); hold on
           %stdshade_acj(squeeze(dvt_data_minus_mean(ch,:,:))',0.2,'g',dvt.times);
           plot(std.times, squeeze(mean(std.data(ch,:,:),2)), 'DisplayName','S'); hold on
           title(chan_name)
           xlabel('Time [ms]');
           ylabel('LFP');
           xline(stim_onset, '--');
           legend('S','D')
           hold off
        end
        %% SAVE ALL
        basename = strcat (participants(i), 'MIs.mat');
        EoI_names = strcat (participants(i), 'EoI.mat');
        cd (basefold)
        c_tag='NC';
        if corrected==1
            c_tag='C';
        end
        cd('MI_Figures') % Save MI figures for the fly

        saveas(h,strcat(char(participants(i)),'_', char(activity_tag), '_',...
            c_tag, num2str(standard_group_number),...
            num2str(deviant_group_number),'_MI_figures_local'),'fig');
        saveas(H,strcat(char(participants(i)),'_', char(activity_tag), '_',...
            c_tag, num2str(standard_group_number),...
            num2str(deviant_group_number),'_ERPs_local'),'fig');
        cd ../
        MI_name = char(strcat(participants(i),char(activity_tag),char(all_con(con)),'_MI_data.mat'));
        save (MI_name,'MI_stat','-mat')

        
        cd('DataEoI') % Save CoI data for the fly
        filename = strcat('EoI_data_', datatype,'.mat');
        if exist(filename,'file')
            load("EoI_data_Drosophila_LFP.mat","EoI")
        end
        EoI.(char(participants(i))).(char(activity_tag)).(char(all_con(con))) = Electorodes;
        EoI.(char(participants(i))).(char(activity_tag)).(char(all_con(con))) = fieldnames(EoI.(char(participants(i))).(char(activity_tag)).(char(all_con(con))));
        save (filename,'EoI','-mat');

        cd ../
        clear Electorodes
        clear MI_stat

    end
end
end




