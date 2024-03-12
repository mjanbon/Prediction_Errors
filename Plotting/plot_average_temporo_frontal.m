%Plots temporo-frontal electrode
function [Temp_front_figure_all] = plot_average_temporo_frontal(basefold, participant, participants, condition, Temporo_front, times, xlimits, ylimits_CoI, ylimit_MI, xticks_CoI, yticks_CoI, yticks_MI, x_labels, y_labels_CoI, y_labels_MI, climits, climits_mask, tempfront_nb)

%%Define matrices for averaged CoI + masks
tempfrontFFiall = [];
tempfrontFFm = zeros(length(times),length(times));
tempfrontFFmr = zeros(length(times),length(times));
tempfrontFFms = zeros(length(times),length(times));
tempfrontFFmi1 = [];
tempfrontFFmi2 = [];
channel_nb = 0;


%% Get data into matrices for plotting

for participanti = 1:length(participants)
   
    tempfrontFFiall = squeeze(cat(3,tempfrontFFiall, Temporo_front.(char(participants(participanti))).data));
    tempfrontFFm = tempfrontFFm + (Temporo_front.(char(participants(participanti))).sigMask / Temporo_front.(char(participants(participanti))).nbchan);
    tempfrontFFmr = tempfrontFFmr + (Temporo_front.(char(participants(participanti))).sigMaskR/ Temporo_front.(char(participants(participanti))).nbchan);
    tempfrontFFms = tempfrontFFms + (Temporo_front.(char(participants(participanti))).sigMaskS/ Temporo_front.(char(participants(participanti))).nbchan);
    tempfrontFFmi1 = squeeze(cat(2,tempfrontFFmi1, Temporo_front.(char(participants(participanti))).MI1));
    tempfrontFFmi2 = squeeze(cat(2,tempfrontFFmi2, Temporo_front.(char(participants(participanti))).MI2));

    participants(participanti)
    Temporo_front.(char(participants(participanti))).nbchan

end

frontFFi = squeeze(mean(tempfrontFFiall,3));
frontFFir = frontFFi;
frontFFis = frontFFi;
frontFFir(frontFFi <= 0 ) = 0;
frontFFis(frontFFi >= 0 ) = 0;


%% PLOT FRONTAL AVERAGE
Temp_front_figure_all = figure(6);
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
stdshade_acj(tempfrontFFmi1',.2,'g',times);
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%MI2
ax(14) = nexttile(4);
stdshade_acj(tempfrontFFmi2',.2,'g',times);  view(-90,90)
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%Combined synergetic/redundant mask
ax(2) = nexttile(6);
imagesc(times,times,tempfrontFFm);
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
imagesc(times,times,tempfrontFFmr);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(6), flipud(bone))
colorbar;

%Synergy mask
ax(7) = nexttile(12);
imagesc(times,times,tempfrontFFms);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(7), flipud(bone))

ax(1) = nexttile(1);
set(gca,'ytick',[], 'yticklabel', [], 'xtick',[], 'xticklabel', [], 'clim', climits);
colorbar; colormap(ax(1), redblue);

filename = char(strcat('Temporo_frontal_average', condition));

cd (strcat(basefold,'Results\Marmo_EcoG\Figures'))
saveas(Temp_front_figure_all,strcat(filename,'.pdf'),'pdf');
saveas(Temp_front_figure_all,strcat(filename,'.fig'),'fig');


end