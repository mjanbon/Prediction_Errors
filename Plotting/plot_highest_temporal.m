%Plots highest temporal electrode
function [Temp_figure, Front_figure] = plot_highest_temporal(basefold, participant, participants, condition, HighAF, times, xlimits, ylimits_CoI, ylimit_MI, xticks_CoI, yticks_CoI, yticks_MI, x_labels, y_labels_CoI, y_labels_MI, climits, climits_mask)

%%Define matrices for averaged CoI + masks
tempFFiall = [];
tempFFm = zeros(length(times),length(times));
tempFFmr = zeros(length(times),length(times));
tempFFms = zeros(length(times),length(times));
tempFFmi1 = [];

frontFFiall = [];
frontFFm = zeros(length(times),length(times));
frontFFmr = zeros(length(times),length(times));
frontFFms = zeros(length(times),length(times));
frontFFmi1 = [];

channel_nb = length(participants);

%% Get data into matrices for plotting

for participanti = 1
    tempFFiall = squeeze(cat(3,tempFFiall, HighAF.(char(participants(participanti))).temp.data));
    tempFFm = tempFFm + logical(HighAF.(char(participants(participanti))).temp.sigMask);
    tempFFmr = tempFFmr + logical(HighAF.(char(participants(participanti))).temp.sigMask > 0);
    tempFFms = tempFFms + logical(HighAF.(char(participants(participanti))).temp.sigMask < 0);
    tempFFmi1 = squeeze(cat(2,tempFFmi1, HighAF.(char(participants(participanti))).temp.MI));

    frontFFiall = squeeze(cat(3,frontFFiall, HighAF.(char(participants(participanti))).front.data));
    frontFFm = frontFFm + logical(HighAF.(char(participants(participanti))).front.sigMask);
    frontFFmr = frontFFmr + logical(HighAF.(char(participants(participanti))).front.sigMask > 0);
    frontFFms = frontFFms + logical(HighAF.(char(participants(participanti))).front.sigMask < 0);
    frontFFmi1 = squeeze(cat(2,frontFFmi1, HighAF.(char(participants(participanti))).front.MI));

end

tempFFi = squeeze(mean(tempFFiall,3));
tempFFir = tempFFi;
tempFFis = tempFFi;
tempFFir(tempFFi <= 0 ) = 0; 
tempFFis(tempFFi >= 0 ) = 0; 

frontFFi = squeeze(mean(frontFFiall,3));
frontFFir = frontFFi;
frontFFis = frontFFi;
frontFFir(frontFFi <= 0 ) = 0; 
frontFFis(frontFFi >= 0 ) = 0; 

frontFFmi1 = cat(2,frontFFmi1,frontFFmi1);
tempFFmi1 = cat(2,tempFFmi1,tempFFmi1);


%% PLOT FRONTAL AVERAGE
Front_figure = figure(6);
tiledlayout(4,3,'TileSpacing','Compact');

%This is the combined synergetic/redundant
ax(1) = nexttile(5);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(times,times,frontFFi, 50,'linecolor','none');
hold on
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(1), redblue)

%MI1
ax(12) = nexttile(2);
stdshade_acj(frontFFmi1',.2,'g',times);
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%MI2
ax(14) = nexttile(4);
stdshade_acj(frontFFmi1',.2,'g',times);  view(-90,90)
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%Combined synergetic/redundant mask
ax(2) = nexttile(6);
imagesc(times,times,frontFFm./channel_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(2), flipud(bone))

%Redundancy CoI
ax(3) = nexttile(8);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(times,times,frontFFir, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(3), redblue)

%synergy COI
ax(4) = nexttile(11);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(times,times,frontFFis, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(4), redblue)

%Redundant mask
ax(6) = nexttile(9);
imagesc(times,times,frontFFmr./channel_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(6), flipud(bone))
colorbar;

%Synergy mask
ax(7) = nexttile(12);
imagesc(times,times,frontFFms./channel_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(7), flipud(bone))

ax(1) = nexttile(1);
set(gca,'ytick',[], 'yticklabel', [], 'xtick',[], 'xticklabel', [], 'clim', clim);
colorbar; colormap(ax(1), redblue);

filename = char(strcat(participant,'highest_average_frontal_', condition));

cd (strcat(basefold,'Results\Marmo_EcoG\Figures'))
saveas(Front_figure,strcat(filename,'.pdf'),'pdf');
saveas(Front_figure,strcat(filename,'.fig'),'fig');

%% PLOT TEMPORAL AVERAGE
Temp_figure = figure(7);
tiledlayout(4,3,'TileSpacing','Compact');

%This is the combined synergetic/redundant
ax(1) = nexttile(5);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(times,times,tempFFi, 50,'linecolor','none');
hold on
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(1), redblue)

%MI1
ax(12) = nexttile(2);
stdshade_acj(tempFFmi1',.2,'g',times);
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%MI2
ax(14) = nexttile(4);
stdshade_acj(tempFFmi1',.2,'g',times);  view(-90,90)
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%Combined synergetic/redundant mask
ax(2) = nexttile(6);
imagesc(times,times,tempFFm./channel_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(2), flipud(bone))

%Redundancy CoI
ax(3) = nexttile(8);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(times,times,tempFFir, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(3), redblue)

%synergy COI
ax(4) = nexttile(11);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(times,times,tempFFis, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(4), redblue)

%Redundant mask
ax(6) = nexttile(9);
imagesc(times,times,tempFFmr./channel_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(6), flipud(bone))
colorbar;

%Synergy mask
ax(7) = nexttile(12);
imagesc(times,times,tempFFms./channel_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(7), flipud(bone))

ax(1) = nexttile(1);
set(gca,'ytick',[], 'yticklabel', [], 'xtick',[], 'xticklabel', [], 'clim', clim);
colorbar; colormap(ax(1), redblue);

filename = char(strcat(participant,'highest_average_temporal_', condition));

cd (strcat(basefold,'Results\Marmo_EcoG\Figures'))
saveas(Temp_figure,strcat(filename,'.pdf'),'pdf');
saveas(Temp_figure,strcat(filename,'.fig'),'fig');

end