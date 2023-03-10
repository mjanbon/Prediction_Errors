%% PLOTS COI/MI WITH PERMUTATION TESTING FOR SPECIFIED PARTICIPANT AND CONDITION
%...in the way we want to have it in the papers 

%Plot JI_XY
%Either within or between CoI

%% GENERAL

% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

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

%PLOTTING PARAMETERS!
cutoff = 16;
xlimits = [700 1150];
ylimits_CoI = [700 1150];
ylimit_MI = [0 0.1];
xticks_CoI = 800:150:110;
yticks_CoI = 800:150:110;
yticks_MI = [0 0.025 0.05 0.075 0.1];
x_labels = 0:150:300;
y_labels_CoI = 0:150:300;
y_labels_MI = [0 0.025 0.05 0.075 0.1];
climits = [-0.01 0.01];
climits_mask = [0 1];

%standard parameters
[basefold, datatype, participantnum, ~, condition, participants, ~, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, ~, start_cut_off, end_cut_off, kperm] = Get_param(0);

%import data to get timing vector
[iEEG1, iEEG2] = impiEEG(participantnum, basefold, datatype, condition, srate, low_cutoff, high_cutoff, filt_order, re_epoch, dev_epochs, std_epochs, epoch_length);

%get proper timing
iEEG1.times (start_cut_off) = [];
iEEG1.times (end_cut_off:end) = [];

%Get data in correct format
[tempFFiall, tempFFm, tempFFmr, tempFFms, tempFFmi1, tempFFmi2, tempFFir, tempFFis, temp_nb, front_nb] = get_plotting_CoI(datatype, char(participants(participantnum)), condition, cutoff, iEEG1.times);

%% Plot Co-I, MI and Masks
a = figure(1);
tiledlayout(4,3,'TileSpacing','Compact');

%This is the combined synergetic/redundant
ax(1) = nexttile(5);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
timing = dvt.times;
contourf(timing,timing,FFi, 50,'linecolor','none');
hold on
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', clim);
colormap(ax(1), redblue)

%MI1
ax(12) = nexttile(2);
stdshade_acj(FFmi1',.2,'g',timing);
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%MI2
ax(14) = nexttile(4);
stdshade_acj(FFmi2',.2,'g',timing);  view(-90,90)
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%Combined synergetic/redundant mask
ax(2) = nexttile(6);
imagesc(timing,timing,FFm./temp_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(2), flipud(bone))

%Redundancy CoI
ax(3) = nexttile(8);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,FFir, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(3), redblue)

%synergy COI
ax(4) = nexttile(11);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,FFis, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(4), redblue)

%Redundant mask
ax(6) = nexttile(9);
imagesc(timing,timing,FFmr./temp_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(6), flipud(bone))
colorbar;

%Synergy mask
ax(7) = nexttile(12);
imagesc(timing,timing,FFms./temp_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(7), flipud(bone))

ax(1) = nexttile(1);
set(gca,'ytick',[], 'yticklabel', [], 'xtick',[], 'xticklabel', [], 'clim', clim);
colorbar; colormap(ax(1), redblue);

saveas(a,'VT_BB_within_local.pdf','pdf')
saveas(a,'VT_BB_within_local.png','png')

