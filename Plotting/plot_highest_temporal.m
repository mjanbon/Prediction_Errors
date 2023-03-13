%Plots highest temporal electrode
function [CoI_figure] = plot_highest_temporal(basefold, participant, participants, condition, HighAF, timing, xlimits, ylimits_CoI, ylimit_MI, xticks_CoI, yticks_CoI, yticks_MI, x_labels, y_labels_CoI, y_labels_MI, climits, climits_mask)

%%Define matrices for averaged CoI + masks
tempFFiall = [];
tempFFm = zeros(length(times),length(times));
tempFFmr = zeros(length(times),length(times));
tempFFms = zeros(length(times),length(times));
tempFFmi1 = [];
tempFFmi2 = [];

frontFFiall = [];
frontFFm = zeros(length(times),length(times));
frontFFmr = zeros(length(times),length(times));
frontFFms = zeros(length(times),length(times));
frontFFmi1 = [];
frontFFmi2 = [];

%% Get data into matrices for plotting

for participanti = 1:length(participants)
    elec_name = char(strcat(temp_elecs(temp_eleci),'_',temp_elecs(temp_eleci)));
    tempFFiall = squeeze(cat(3,tempFFiall,CoI.(partname).(condition).data.(elec_name)));
    tempFFm = tempFFm + logical(CoI.(partname).(condition).sigMask.(elec_name));
    temp_all_mask = squeeze(cat(3,temp_all_mask,CoI.(partname).(condition).sigMask.(elec_name)));
    tempFFmr = tempFFmr + logical(CoI.(partname).(condition).sigMask.(elec_name)>0);
    tempFFms = tempFFms + logical(CoI.(partname).(condition).sigMask.(elec_name)<0);
    tempFFmi1 = squeeze(cat(2,tempFFmi1, CoI.(partname).(condition).MI1.(elec_name)));
    tempFFmi2 = squeeze(cat(2,tempFFmi2, CoI.(partname).(condition).MI2.(elec_name)));
end


%% Create figure and plot
CoI_figure = figure(1);
tiledlayout(4,3,'TileSpacing','Compact');

%This is the combined synergetic/redundant
ax(1) = nexttile(5);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,frontintFFi, 50,'linecolor','none');
hold on
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(1), redblue)

%MI1
ax(12) = nexttile(2);
stdshade_acj(frontintFFmi1',.2,'g',timing);
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%MI2
ax(14) = nexttile(4);
stdshade_acj(frontintFFmi2',.2,'g',timing);  view(-90,90)
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);

%Combined synergetic/redundant mask
ax(2) = nexttile(6);
imagesc(timing,timing,frontintFFm./frontint_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(2), flipud(bone))

%Redundancy CoI
ax(3) = nexttile(8);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,frontintFFir, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(3), redblue)

%synergy COI
ax(4) = nexttile(11);
set (gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,frontintFFis, 50,'linecolor','none');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits);
colormap(ax(4), redblue)

%Redundant mask
ax(6) = nexttile(9);
imagesc(timing,timing,frontintFFmr./frontint_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(6), flipud(bone))
colorbar;

%Synergy mask
ax(7) = nexttile(12);
imagesc(timing,timing,frontintFFms./frontint_nb);
set(gca,'YDir','normal');
xlim(xlimits);  ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, 'xtick',xticks_CoI, 'xticklabel', x_labels, 'clim', climits_mask);
colormap(ax(7), flipud(bone))

ax(1) = nexttile(1);
set(gca,'ytick',[], 'yticklabel', [], 'xtick',[], 'xticklabel', [], 'clim', clim);
colorbar; colormap(ax(1), redblue);

filename = char(strcat(participant_name,'_temporo_temporal_', condition));

cd (strcat(basefold,'Results\Marmo_EcoG\Figures'))
saveas(CoI_figure,strcat(filename,'.pdf'),'pdf');
saveas(CoI_figure,strcat(filename,'.fig'),'fig');