%% ESTIMATES MI AND GETS ELECTRODES OF INTEREST
% Estimates MI via gaussian copula estimation (Ince et al., 2016), performs
% parametric permutation testing and saves electrodes with significant MI
% into a mat-file.


function main_diff_MI_ERP(task_id)
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

[basefold, datatype, all_con, condition, subject,participants, EoI,...
    srate, activity_tag, ~, ~, corrected, stim_onset, baseline,...
    start_cut_off, end_cut_off, kperm] = Max_get_param(USING_HPC, 0);
kperm=100
activity_tags = {'wake', 'sleep'}
%% Load the group_numbers manually for sleep and wake
deviant_group_number_wake = [7, 8]; %Deviant group in groupHyper for decomposed sleep/wake (wake)
standard_group_number_wake = [11, 12]; %Carrier group in groupHyper for decomposed sleep/wake (wake)

deviant_group_number_sleep = [3, 4]; %Deviant group in groupHyper for decomposed sleep/wake (all sleep)
standard_group_number_sleep = [9, 10]; %Carrier group in groupHyper for decomposed sleep/wake (all sleep)

%Go through all specified participants and conditions
for i = 1: length(participants)
    for con = 1:length(all_con)
        %% IMPORT DATA
        participants(i)
        all_con(con)
        % Import data and check for same amount of trials & channels
        % [dvt, std] = impiEEG(i, basefold, datatype, all_con(con), srate, low_cutoff, high_cutoff, filt_order,re_epoch, dev_epochs, std_epochs, epoch_length);
        overVar_file = strcat(basefold, datatype, '/', participants(i), '_', all_con(con), '.mat')
        if exist(char(overVar_file)) == 0
            continue
        end

        [dvt_w, std_w] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_wake, standard_group_number_wake,...
            corrected, srate); % Load the wake deviants and standards
        
        [dvt_s, std_s] = load_trials_from_group_hyper(char(overVar_file), deviant_group_number_sleep, standard_group_number_sleep,...
            corrected, srate); % Load the wake deviants and standards
        % 1. Match number of trials for sleep and wake (deviant and standard)
        n_trials_dvt = min([dvt_w.trials, dvt_s.trials]);
        n_trials_std = min([std_w.trials, std_s.trials]);

        dvt_w.data = dvt_w.data(:,:,1:n_trials_dvt); size(dvt_w.data)
        dvt_s.data = dvt_s.data(:,:,1:n_trials_dvt); size(dvt_s.data)
        std_w.data = std_w.data(:,:,1:n_trials_std); size(dvt_w.data)
        std_s.data = std_s.data(:,:,1:n_trials_std); size(dvt_s.data)

        % 2. Baseline normalization (for deviant and standard, sleep and wake)
        % (Assume 'baseline' is a vector of baseline sample indices)
        for cond = {'dvt_w','dvt_s','std_w','std_s'}
            dat = eval(cond{1});
            data = permute(dat.data, [1 3 2]); % channels x trials x samples
            mean_data = squeeze(mean(data(:,:,baseline),3));
            dat.data = data - repmat(mean_data, [1 1 size(dat.data,2)]);
            dat.data = permute(dat.data, [1 3 2]); % back to channels x samples x trials
            eval([cond{1} ' = dat;']);
        end
        
        % Prepare save directory
        save_dir = fullfile(basefold, 'MI_figures', 'wake_sleep_permutations');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end

        % Create a single figure for all channels
        figure('Position', [100 100 1800 900]);
        tiledlayout(3, 5, 'TileSpacing', 'compact', 'Padding', 'compact'); % adjust for your number of channels

        for ch = 1:dvt_w.nbchan
            % Wake and sleep data: [samples x trials]
            data_wake = [squeeze(dvt_w.data(ch,:,:))'; squeeze(std_w.data(ch,:,:))'];
            data_sleep = [squeeze(dvt_s.data(ch,:,:))'; squeeze(std_s.data(ch,:,:))'];
            stim = [ones(n_trials_dvt,1); zeros(n_trials_std,1)];

            % Ensure correct orientation: [samples x trials]
            if size(data_wake,2) ~= length(stim)
                data_wake = data_wake';
            end
            if size(data_sleep,2) ~= length(stim)
                data_sleep = data_sleep';
            end

            nPerm = kperm;
            [MI_diff, MI_wake, MI_sleep, sigMask] = cnm_MI_stimtime_fwer_diff(data_wake, data_sleep, stim, nPerm);

            % Plot in tiled layout
            nexttile;
            timing = dvt_w.times;
            plot(timing, MI_wake, 'r', 'DisplayName', 'Wake'); hold on;
            plot(timing, MI_sleep, 'b', 'DisplayName', 'Sleep');
            plot(timing, MI_diff, 'k', 'DisplayName', 'Wake - Sleep');
            sig_idx = find(sigMask ~= 0);
            if ~isempty(sig_idx)
                scatter(timing(sig_idx), MI_diff(sig_idx), 40, 'm', 'filled', 'DisplayName', 'Significant');
            end
            title(['Ch ' num2str(ch)]);
            xlabel('Time [ms]');
            ylabel('MI (bits)');
            if ch == 1
                legend;
            end
            hold off;
        end

        sgtitle(['MI difference (Wake - Sleep) for all channels']);

        % Save the figure
        saveas(gcf, fullfile(save_dir, [char(participants(i)) '_' char(all_con(con)) '_MI_diff_all_channels.fig']));
        saveas(gcf, fullfile(save_dir, [char(participants(i)) '_' char(all_con(con)) '_MI_diff_all_channels.png']));
        print(gcf, fullfile(save_dir, [char(participants(i)) '_' char(all_con(con)) '_MI_diff_all_channels.svg']), '-dsvg');
        close(gcf);
    end

            % Optionally, save results
            % save(fullfile(save_dir, sprintf('MI_diff_ch%d.mat', ch)), 'MI_diff', 'MI_wake', 'MI_sleep', 'sigMask', 'timing');
        % %% BASELINE NORMALISATION FOR WAKE AND SLEEP
        % data_dvt_w = permute(dvt_w.data,[1 3 2]); mean_data_dvt_w = squeeze(mean(data_dvt_w(:,:,baseline),3));
        % dvt_w.data = data_dvt - repmat(mean_data_dvt_w,[1 1 size(dvt.data_w,2)]);
        % bb1_dev_bl = permute(dvt_w.data,[1 3 2]);

        % data_std = permute(std.data,[1 3 2]); mean_data_std = squeeze(mean(data_std(:,:,baseline),3));
        % std.data = data_std - repmat(mean_data_std,[1 1 size(std.data,2)]);
        % % std.data(:,:,start_cut_off) = [];
        % % std.data(:,:,end_cut_off:end) = [];
        % bb1_std_bl = permute(std.data,[1 3 2]);

        % bb_dev = permute(bb1_dev_bl,[3 1 2]);
        % bb_std = permute(bb1_std_bl,[3 1 2]);

        % %% CALCULATE MUTUAL INFORMATION & PERFORM PERMUTATION TESTING USING GCMI

        % elec_of_I = {};
        % h = figure(con);
        % for ch = 1 : dvt.nbchan
        %     chnum = strcat('ch_',string(ch));
        %     dat.class1.BB = squeeze(bb_dev(:, ch, :));
        %     dat.class2.BB = squeeze(bb_std(:, ch, :));

        %     %Calculate MI between the specified channels
        %     [MI, sigMask] = cnm_MI_stimtime([dat.class1.BB; dat.class2.BB],[zeros(1, size(bb_dev,1)), ones(1, size(bb_dev,1))]', kperm);

        %     %save
        %     chan_name = strcat('E',num2str(ch));
        %     MI_stat.(participants{i}).(all_con{con}).MI.(chan_name) = MI;
        %     MI_stat.(participants{i}).(all_con{con}).electrode = chan_name;
        %     MI_stat.(participants{i}).(all_con{con}).sigMask.(chan_name)= sigMask;

        %     %% PLOT
        %     %get correct timing
        %     timing = dvt.times;
        %     % timing(start_cut_off) = [];
        %     % timing(end_cut_off:end) = [];
        %     % timing = timing';

        %     %plot MI
        %     nexttile(ch)
        %     plot(timing, MI_stat.(participants{i}).(all_con{con}).MI.(chan_name));
        %     title(chan_name)
        %     hold on;
        %     xlabel('Time [ms]');
        %     yLimits = [-0.001 0.2];
        %     %         xlim([0.7 1.150])
        %     pos_sigbar = yLimits(2) - (0.07 * range(yLimits));
        %     xline(stim_onset, '--');


        %     % Draw stats
        %     stat= logical(MI_stat.(participants{i}).(all_con{con}).sigMask.(chan_name));
        %     ylim(yLimits);
        %     ylabel('Mutual information (bits)');

        %     EoI_list = false;
        %     for tIdx = 1:length(timing)-1
        %         tIdx2 = timing(tIdx);
        %         tIdx2_1 = timing(tIdx+1);
        %         if stat(tIdx) > 0
        %             plot([tIdx2, tIdx2_1], [pos_sigbar, pos_sigbar], 'LineWidth', 3, 'Color', 'm');
        %             EoI_list = true;
        %         end
        %     end

        %     %List electrodes of interest
        %     if EoI_list == true
        %         chan = (MI_stat.(participants{i}).(all_con{con}).electrode); %Not electrode(i) anymore
        %         chan = replace(chan,'-','_');
        %         Electorodes.(chan_name) = chan_name;
        %     end
        %     hold off
        % end
        % %
        % %% PLOT ERPS

        % % dvt_data_minus_mean = permute(dvt_data_minus_mean,[2 3 1]);
        % % std_data_minus_mean = permute(std_data_minus_mean,[2 3 1]);
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
        % %% SAVE ALL
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
        % cd MI_Data
        % MI_name = char(strcat(participants(i),char(activity_tag),char(all_con(con)),'_MI_data.mat'));
        % save (MI_name,'MI_stat','-mat')
        % cd ../

        
        % cd('DataEoI') % Save CoI data for the fly
        % filename = strcat('EoI_data_', datatype,'.mat');
        % if exist(filename,'file')
        %     load("EoI_data_Drosophila_LFP.mat","EoI")
        % end
        % EoI.(char(participants(i))).(char(activity_tag)).(char(all_con(con))) = Electorodes;
        % EoI.(char(participants(i))).(char(activity_tag)).(char(all_con(con))) = fieldnames(EoI.(char(participants(i))).(char(activity_tag)).(char(all_con(con))));
        % save (filename,'EoI','-mat');

        % cd ../
        % clear Electorodes
        % clear MI_stat
        % close all


end





