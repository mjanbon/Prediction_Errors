%Plots tempo-frontal electrodes 
%saves figures as Pdf and fig
function [CoI_figure] = plot_tempo_frontal(basefold,participant_name, condition, tempfrontFFi, tempfrontFFm, tempfrontFFmr, tempfrontFFms, tempfrontFFmi1, tempfrontFFmi2, tempfrontFFir, tempfrontFFis, tempfront_nb, timing, ...
                                        xlimits,ylimits_CoI,ylimit_MI,xticks_CoI,yticks_CoI,yticks_MI,x_labels,y_labels_CoI, y_labels_MI,climits,climits_mask)

%% Create figure and plot
CoI_figure = figure(1);
tiledlayout(4,3,'TileSpacing','Compact');

%This is the combined synergetic/redundant
ax(1) = nexttile(5);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempfrontFFi, 50,'linecolor','none');
hold on
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(1), redblue)

%MI1
ax(12) = nexttile(2);
stdshade_acj(tempfrontFFmi1',.2,'g',timing);
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%MI2
ax(14) = nexttile(4);
stdshade_acj(tempfrontFFmi2',.2,'g',timing);  view(-90,90)
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%Combined synergetic/redundant mask
ax(2) = nexttile(6);
imagesc(timing,timing,tempfrontFFm./tempfront_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(2), flipud(bone))

%Redundancy CoI
ax(3) = nexttile(8);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempfrontFFir, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(3), redblue)

%synergy COI
ax(4) = nexttile(11);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempfrontFFis, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(4), redblue)

%Redundant mask
ax(6) = nexttile(9);
imagesc(timing,timing,tempfrontFFmr./tempfront_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(6), flipud(bone))
colorbar;

%Synergy mask
ax(7) = nexttile(12);
imagesc(timing,timing,tempfrontFFms./tempfront_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(7), flipud(bone))

ax(1) = nexttile(1);
set(gca,'ytick',[], 'yticklabel', [], 'xtick',[], 'xticklabel', [], 'clim', clim);
colorbar; colormap(ax(1), redblue);

filename = char(strcat(participant_name,'_temporo_frontal_', condition));

cd (strcat(basefold,'Results\Marmo_EcoG\Figures'))
saveas(CoI_figure,strcat(filename,'.pdf'),'pdf');
saveas(CoI_figure,strcat(filename,'.fig'),'fig');