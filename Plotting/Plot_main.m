%% PLOTS COI/MI WITH PERMUTATION TESTING FOR SPECIFIED PARTICIPANT AND CONDITION
%...in the way we want to have it in the papers :)

%% GENERAL
clear all
% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperSize=[29.7 29.7];
fig.PaperPosition = [1 1 20 28];
fig.PaperType = 'a4' ;
fig.PaperOrientation = 'portrait' ;

% general graphics, this will apply to any figure you open (groot is the default figure object).
% I have this in my startup.m file, so I don't have to retype these things whenever plotting a new fig.
set(groot, ...
    'DefaultFigureColorMap', linspecer, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 0.5, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 8, ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 8, ...
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.025]);

% type cbrewer without input args to see all possible sets of colormaps
colors = cbrewer2('qual', 'Set1', 8);

%% GETS ALL THE PARAMETERS FOR PARTICIPANT AND DATA

%PARTICIPANT PARAMETERS (CHANGE IN THE Get_param-function)
[basefold, datatype, participanti, ~, condition, participants, ~, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, ~, start_cut_off, end_cut_off, kperm] = Get_param(0);

EoI_file = strcat('EoI_data_',datatype,'.mat');
load(EoI_file)

conditions = {'attention' 'counting'}
for participanti = 11
for conditioni = 1:2
    condition = char(conditions(conditioni));

%     if strcmp(condition,'counting') == 1 && strcmp(char(participants(participanti)),'EOW')
%         continue
%     end

cut1 = EoI.(char(participants(participanti))).cutoff_temp.attention +1;
cut2 = EoI.(char(participants(participanti))).cutoff_temp.counting +1;

% temporal_attention = EoI.AML.attention(1:cut1);
% temporal_counting = EoI.AML.counting(1:cut1);
% 
% frontal_attention = EoI.AML.attention(cut1:end);
% frontal_counting = EoI.AML.counting(cut1:end);

%PLOTTING PARAMETERS!
cut = [cut1, cut2];
xlimits = [-100 600];
ylimits_CoI = [-100 600];
ylimit_MI = [-0.005 0.1];
xticks_CoI = 0:150:600
yticks_CoI = 0:150:600;
yticks_MI = [0:0.05:0.1]
x_labels = 0:150:600;
y_labels_CoI = 0:150:600;
y_labels_MI = [0:0.05:0.1];
climits = [-0.01 0.01];
climits_mask = [0 1];

%% LOAD DATA

%import data to get timing vector

timing = -99:2:600;
% timing = -299:2:800;

%Get data in correct format

% %FOR PLOTTING AVERAGE OF ALL
% [tempFFiall, tempFFi, temp_all_mask, tempFFm, tempFFmr, tempFFms,tempFFmi1, tempFFmi2,tempFFir,tempFFis, temp_nb, front_nb, tempfront_nb, tempint_nb, frontint_nb, ...
%           frontFFiall, frontFFi, front_all_mask, frontFFm, frontFFmr, frontFFms, frontFFmi1, frontFFmi2, frontFFir, frontFFis,...
%           tempfrontFFiall, tempfrontFFi, tempfrontFFm, tempfrontFFmr, tempfrontFFms, tempfrontFFmi1, tempfrontFFmi2, tempfrontFFir, tempfrontFFis,...
%           tempintFFiall, tempintFFi, tempintFFm, tempintFFmr, tempintFFms, tempintFFmi1, tempintFFmi2, tempintFFir, tempintFFis,...
%           frontintFFiall, frontintFFi, frontintFFm, frontintFFmr, frontintFFms, frontintFFmi1, frontintFFmi2, frontintFFir, frontintFFis, temp_elecs, front_elecs] = get_plotting_CoI(basefold,datatype, participanti, char(participants(participanti)), condition, cut, timing);

%%FOR PLOTTING individual
[tempFFiall, tempFFi, temp_all_mask, tempFFm, tempFFmr, tempFFms,tempFFmi1, tempFFmi2,tempFFir,tempFFis, temp_nb, front_nb, tempfront_nb, tempint_nb, frontint_nb, ...
          frontFFiall, frontFFi, front_all_mask, frontFFm, frontFFmr, frontFFms, frontFFmi1, frontFFmi2, frontFFir, frontFFis,...
          tempfrontFFiall, tempfrontFFi, tempfrontFFm, tempfrontFFmr, tempfrontFFms, tempfrontFFmi1, tempfrontFFmi2, tempfrontFFir, tempfrontFFis,...
          tempintFFiall, tempintFFi, tempintFFm, tempintFFmr, tempintFFms, tempintFFmi1, tempintFFmi2, tempintFFir, tempintFFis,...
          frontintFFiall, frontintFFi, frontintFFm, frontintFFmr, frontintFFms, frontintFFmi1, frontintFFmi2, frontintFFir, frontintFFis, temp_elecs, front_elecs] = get_plotting_CoI_ALL(basefold,datatype, participanti, char(participants(participanti)), condition, cut, timing);



%% % PLOT Averages of all electrodes ( temporal/frontal/fronto-temporal/temporo-temporal/fronto-frontal)
% % 
%Plot temporal figures
if isempty(tempFFi) == 0
[temp_figure] = plot_temporal(basefold,participants(participanti), condition, tempFFi, tempFFm, tempFFmr, tempFFms, tempFFmi1, tempFFmi2, tempFFir, tempFFis, temp_nb, timing,...
    xlimits,ylimits_CoI,ylimit_MI,xticks_CoI,yticks_CoI,yticks_MI,x_labels,y_labels_CoI, y_labels_MI,climits,climits_mask);
end 
% 
% if isempty(frontFFi) == 0
% % Plot frontal figures
% [front_figure] = plot_frontal(basefold,participants(participanti), condition, frontFFi, frontFFm, frontFFmr, frontFFms, frontFFmi1, frontFFmi2, frontFFir, frontFFis, front_nb, timing,...
%     xlimits,ylimits_CoI,ylimit_MI,xticks_CoI,yticks_CoI,yticks_MI,x_labels,y_labels_CoI, y_labels_MI,climits,climits_mask);
% end
% 
% if isempty(tempfrontFFi) == 0
% % Plot temporo-frontal figures
% [tempo_front_figure] = plot_tempo_frontal(basefold,participants(participanti), condition, tempfrontFFi, tempfrontFFm, tempfrontFFmr, tempfrontFFms, tempFFmi1, frontFFmi1, tempfrontFFir, tempfrontFFis, tempfront_nb, timing,...
%     xlimits,ylimits_CoI,ylimit_MI,xticks_CoI,yticks_CoI,yticks_MI,x_labels,y_labels_CoI, y_labels_MI,climits,climits_mask);
% end
% 
% if isempty(tempintFFi) == 0
% % Plot temporo-temporal figures
% [temporo_temporal_figure] = plot_tempo_temporal(basefold,participants(participanti), condition, tempintFFi, tempintFFm, tempintFFmr, tempintFFms, tempFFmi1, tempFFmi2, tempintFFir, tempintFFis, tempint_nb, timing,...
%     xlimits,ylimits_CoI,ylimit_MI,xticks_CoI,yticks_CoI,yticks_MI,x_labels,y_labels_CoI, y_labels_MI,climits,climits_mask);
% end
% 
% if isempty(frontintFFi) == 0
% % Plot temporo-temporal figures
% [frontal_frontal_figure] = plot_fronto_frontal(basefold,participants(participanti), condition, frontintFFi, frontintFFm, frontintFFmr, frontintFFms, frontFFmi1, frontFFmi2, frontintFFir, frontintFFis, frontint_nb, timing,...
%     xlimits,ylimits_CoI,ylimit_MI,xticks_CoI,yticks_CoI,yticks_MI,x_labels, y_labels_CoI, y_labels_MI,climits, climits_mask);
% end


% Plot the electrodes with the highest MI
% find the electrodes with the highest MI for all participants and load data into struct
% for participanti = 1:length(participants)
% 
% [tempFFiall, tempFFi, temp_all_mask, tempFFm, tempFFmr, tempFFms,tempFFmi1, tempFFmi2,tempFFir,tempFFis, temp_nb, front_nb, tempfront_nb, tempint_nb, frontint_nb, ...
%           frontFFiall, frontFFi, front_all_mask, frontFFm, frontFFmr, frontFFms, frontFFmi1, frontFFmi2, frontFFir, frontFFis,...
%           tempfrontFFiall, tempfrontFFi, tempfrontFFm, tempfrontFFmr, tempfrontFFms, tempfrontFFmi1, tempfrontFFmi2, tempfrontFFir, tempfrontFFis,...
%           tempintFFiall, tempintFFi, tempintFFm, tempintFFmr, tempintFFms, tempintFFmi1, tempintFFmi2, tempintFFir, tempintFFis,...
%           frontintFFiall, frontintFFi, frontintFFm, frontintFFmr, frontintFFms, frontintFFmi1, frontintFFmi2, frontintFFir, frontintFFis, temp_elecs, front_elecs] = get_plotting_CoI(basefold,datatype, participanti, char(participants(participanti)), condition, cutoff, timing);
% 
% Temporo_front.(char(participants(participanti))).data = tempfrontFFi;
% Temporo_front.(char(participants(participanti))).sigMask = tempfrontFFm;
% Temporo_front.(char(participants(participanti))).sigMaskR = tempfrontFFmr;
% Temporo_front.(char(participants(participanti))).sigMaskS = tempfrontFFms;
% Temporo_front.(char(participants(participanti))).MI1 = tempFFmi1;
% Temporo_front.(char(participants(participanti))).MI2 = frontFFmi1;
% Temporo_front.(char(participants(participanti))).nbchan = tempfront_nb;
% 
% [Highest, MI] = Get_highest_MI(basefold,datatype,condition, participanti, char(participants(participanti)), cut, timing,participants);
% HighAF.(char(participants(participanti))) = Highest;

% timing = -199:2:600;
% Plot highest temporal electrodes for all participants (+ individual plots)
% [temp_high_fig, front_high_fig, temp_high_all, front_high_all]  = plot_highest(basefold, participants(participanti), participants, condition, HighAF, timing, xlimits, ylimits_CoI, ylimit_MI, xticks_CoI, yticks_CoI, yticks_MI, x_labels, y_labels_CoI, y_labels_MI, climits, climits_mask, Highest, participanti);

% if isempty(tempFFmi1) == 0
% temp_MI = figure (8);
% stdshade_acj(tempFFmi1',.2,'g',timing);
% title('temporal MI')
% xlim(xlimits);  ylim(ylimit_MI);
% set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);
% end
% 
% if isempty(frontFFmi1) == 0
% front_MI = figure(9);
% stdshade_acj(frontFFmi1',.2,'g',timing);
% title('Frontal MI')
% xlim(xlimits);  ylim(ylimit_MI);
% set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);
% end

% [front, temp] = Highest_MI(MI,condition,participants(participanti));

part = char(participants(participanti));
cond = condition;
coi= 'CoI_figures';
cd ('D:\iEEG\Data\All_data\Human_attention\')
status = mkdir(coi);
cd ('D:\iEEG\Data\All_data\Human_attention\CoI_figures\')
status = mkdir(part);
cd (strcat('D:\iEEG\Data\All_data\Human_attention\CoI_figures\',part,'\'));
status = mkdir(cond);
cd (strcat('D:\iEEG\Data\All_data\Human_attention\CoI_figures\',part,'\',cond,'\'))
            
            if isempty(tempFFi) == 0
            name = char(strcat(participants(participanti),'_',condition,'_temp_within_bb'));
            saveas(temp_figure, name ,'pdf')
            saveas(temp_figure, name ,'png')
            saveas(temp_figure, name ,'fig')

            clear temp_figure
            end
% 
%             if isempty(frontFFi) == 0
%             name = char(strcat(participants(participanti),'_',condition,'_front_within_bb'));
%             saveas(front_figure, name ,'pdf')
%             saveas(front_figure, name ,'jpg')
%             saveas(front_figure, name ,'fig')
% 
%             clear front_figure
%             end
% 
%             if isempty(tempfrontFFi) == 0
%             name = char(strcat(participants(participanti),'_',condition,'_temporo-frontal_bb'));
%             saveas(tempo_front_figure, name ,'pdf')
%             saveas(tempo_front_figure, name ,'jpg')
%             saveas(tempo_front_figure, name ,'fig')
% 
%             clear tempo_front_figure
%             end
%             
%             if isempty(tempintFFi) == 0
%             name = char(strcat(participants(participanti),'_',condition,'_temporal-between_bb'));
%             saveas(temporo_temporal_figure, name ,'pdf')
%             saveas(temporo_temporal_figure, name ,'jpg')
%             saveas(temporo_temporal_figure, name ,'fig')
%             end
%             
%             if isempty(frontintFFi) == 0
%             name = char(strcat(participants(participanti),'_',condition,'_frontal-between_bb'));
%             saveas(frontal_frontal_figure, name ,'pdf')
%             saveas(frontal_frontal_figure, name ,'jpg')
%             saveas(frontal_frontal_figure, name ,'fig')
% 
%             clear frontal_frontal_figure
%             end
% 
%             if isempty(tempFFmi1) == 0
%             name = char(strcat(participants(participanti),'_',condition,'_temp_MI_bb'));
%             saveas(temp_MI, name ,'pdf')
%             saveas(temp_MI, name ,'jpg')
%             saveas(temp_MI, name ,'fig')
% 
%             clear temp_MI
%             end
%             
%             if isempty(frontFFmi1) == 0
%             name = char(strcat(participants(participanti),'_',condition,'_front_MI_bb'));
%             saveas(front_MI, name ,'pdf')
%             saveas(front_MI, name ,'jpg')
%             saveas(front_MI, name ,'fig')
%             clear front_MI
%             end
% 
%             if isempty(tempFFmi1) == 0
%             name = char(strcat(participants(participanti),'_',condition,'_highest_temp_bb'));
%             saveas(temp_high_fig, name ,'pdf')
%             saveas(temp_high_fig, name ,'jpg')
%             saveas(temp_high_fig, name ,'fig')
% 
%             clear temp_high_fig
%             end
% 
%             if isempty(frontFFmi1) == 0
%             name = char(strcat(participants(participanti),'_',condition,'_highest_front_bb'));
%             saveas(front_high_fig, name ,'pdf')
%             saveas(front_high_fig, name ,'jpg')
%             saveas(front_high_fig, name ,'fig')
% 
%             clear front_high_fig
%             end


%             name = char(strcat(participants(participanti),'_',condition,'_comparison_front'));
%             saveas(front, name ,'pdf')
%             saveas(front, name ,'jpg')
%             saveas(front, name ,'fig')
% 
%             clear front
% 
%             name = char(strcat(participants(participanti),'_',condition,'_comparison_temp'));
%             saveas(temp, name ,'pdf')
%             saveas(temp, name ,'jpg')
%             saveas(temp, name ,'fig')
% 
%             clear temp
% 
%             name = 'Highest_all_temporal';
%             saveas(temp_high_all, name ,'pdf')
%             saveas(temp_high_all, name ,'jpg')
%             saveas(temp_high_all, name ,'fig')
% 
%             name = 'Highest_all_frontal';
%             saveas(front_high_all, name ,'pdf')
%             saveas(front_high_all, name ,'jpg')
%             saveas(front_high_all, name ,'fig')
end
end















