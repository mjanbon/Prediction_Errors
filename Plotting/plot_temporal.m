%Plots temporal electrodes 
%saves figures as Pdf and fig
function [CoI_figure] = plot_temporal(basefold,participant_name, condition, tempFFi, tempFFm, tempFFmr, tempFFms, tempFFmi1, tempFFmi2, tempFFir, tempFFis, temp_nb, timing, ...
                                        xlimits,ylimits_CoI,ylimit_MI,xticks_CoI,yticks_CoI,yticks_MI,x_labels,y_labels_CoI, y_labels_MI,climits,climits_mask)

%% Create figure and plot
CoI_figure = figure(1);
tiledlayout(4,3,'TileSpacing','Compact');

%This is the combined synergetic/redundant
ax(1) = nexttile(5);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempFFi, 50,'linecolor','none');
hold on
xlim(xlimits);  ylim(ylimits_CoI); clim(climits)
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels);
colormap(ax(1), redblue)

lims = clim
%MI1
ax(12) = nexttile(2);
stdshade_acj(tempFFmi1',.2,'g', timing);
% title(strcat(condition,' temporal within all'))
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%MI2
ax(14) = nexttile(4);
stdshade_acj(tempFFmi2',.2,'g',timing);  view(-90,90)
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%Combined synergetic/redundant mask
ax(2) = nexttile(6);
imagesc(timing,timing,tempFFm./temp_nb);
title(num2str(temp_nb));
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI); clim(climits)
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(2), flipud(bone))

%Redundancy CoI
ax(3) = nexttile(8);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempFFir, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI); clim(climits)
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels);
colormap(ax(3), redblue)

%synergy COI
ax(4) = nexttile(11);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempFFis, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels);
colormap(ax(4), redblue)

%Redundant mask
ax(6) = nexttile(9);
imagesc(timing,timing,tempFFmr./temp_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(6), flipud(bone))
colorbar;

%Synergy mask
ax(7) = nexttile(12);
imagesc(timing,timing,tempFFms./temp_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(7), flipud(bone))

ax(1) = nexttile(1);
set(gca,'ytick',[], 'yticklabel', [], 'xtick',[], 'xticklabel', [], 'clim', lims);
colorbar; colormap(ax(1), redblue);
hold off

filename = char(strcat(participant_name,'_temporal_within_', condition));
% 
% cd (strcat(basefold,'Results\Marmo_EcoG\Figures'))
% saveas(CoI_figure,strcat(filename,'.pdf'),'pdf');
% saveas(CoI_figure,strcat(filename,'.fig'),'fig');
