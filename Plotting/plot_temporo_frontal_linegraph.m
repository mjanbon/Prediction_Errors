%Plots highest temporal electrode
function [Temp_front_figure_line] = plot_temporo_frontal_linegraph(basefold, participant, participants, condition, Temporo_front, times, xlimits, ylimits_CoI, ylimit_MI, xticks_CoI, yticks_CoI, yticks_MI, x_labels, y_labels_CoI, y_labels_MI, climits, climits_mask

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
    tempfrontFFm = tempfrontFFm + Temporo_front.(char(participants(participanti))).sigMask;
    tempfrontFFmr = tempfrontFFmr + Temporo_front.(char(participants(participanti))).sigMaskR;
    tempfrontFFms = tempfrontFFms + Temporo_front.(char(participants(participanti))).sigMaskS;
    tempfrontFFmi1 = squeeze(cat(2,tempfrontFFmi1, Temporo_front.(char(participants(participanti))).MI1));
    tempfrontFFmi2 = squeeze(cat(2,tempfrontFFmi2, Temporo_front.(char(participants(participanti))).MI2));
    channel_nb = channel_nb + Temporo_front.(char(participants(participanti))).nbchan;

end

tempfrontFFi = squeeze(mean(tempfrontFFiall,3));
tempfrontFFir = tempfrontFFi;
tempfrontFFis = tempfrontFFi;
tempfrontFFir(tempfrontFFi <= 0 ) = 0;
tempfrontFFis(tempfrontFFi >= 0 ) = 0;


%% PLOT LINEGRAPH

Temp_front_figure_line = figure (1);
tiledlayout(1,3,'TileSpacing','Compact');

nexttile(1)
stdshade_acj(tempfrontFFi,.2,'g',times);
set(gca,'xlim',[700 1150],'ylim',[-0.0025 0.0025],'ytick',[-0.025 0 0.025], 'yticklabel', [ ], 'xtick',[0 150 300], 'xticklabel', [0 150 300]);

% nexttile(2)
% stdshade_acj(tempfrontFFir,.2,'g',times);
% 
% nexttile(2)
% stdshade_acj(tempfrontFFis,.2,'g',times);


end