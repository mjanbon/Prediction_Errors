load('NrXY_MI_data.mat')

MI_temp_nr = MI_stat.Nr.XY.MI.E24;
Sig_temp_nr = MI_stat.Nr.XY.sigMask.E24;

MI_front_nr = MI_stat.Nr.XY.MI.E27;
Sig_front_nr = MI_stat.Nr.XY.sigMask.E27;
% 
load('JiXY_MI_data.mat')

MI_temp_ji = MI_stat.Ji.XY.MI.E5;
Sig_temp_ji = MI_stat.Ji.XY.sigMask.E5;

MI_front_ji = MI_stat.Ji.XY.MI.E18;
Sig_front_ji = MI_stat.Ji.XY.sigMask.E18;

timing = 700:2:1148;

xlimits = [700 1150];
ylimits_CoI = [700 1150];
ylimit_MI = [-0.005 0.1];
xticks_CoI = 800:150:1100;
yticks_CoI = 800:150:1100;
yticks_MI = [0 0.05];
x_labels = 0:150:300;
y_labels_CoI = 0:150:300;
y_labels_MI = [0 0.05];

figure (1)
%NR_temporal
plot(timing, MI_temp_nr );
title('Nr temporal')
hold on;
xlabel('Time [ms]');
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);
pos_sigbar = ylimit_MI(2) - (0.07 * range(yLimits));
stat=logical(Sig_temp_nr);
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
hold off

%NR_temporal
figure (2)
plot(timing, MI_front_nr );
title('Nr frontal')
hold on;
xlabel('Time [ms]');
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);
pos_sigbar = ylimit_MI(2) - (0.07 * range(yLimits));
stat=logical(Sig_front_nr);
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

figure (3)
%NR_temporal
plot(timing, MI_temp_ji );
title('Ji temporal')
hold on;
xlabel('Time [ms]');
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);
pos_sigbar = ylimit_MI(2) - (0.07 * range(yLimits));
stat=logical(Sig_temp_ji);
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
hold off

%NR_temporal
figure (4)
plot(timing, MI_front_ji );
title('Ji frontal')
hold on;
xlabel('Time [ms]');
xlim(xlimits);  ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);
pos_sigbar = ylimit_MI(2) - (0.07 * range(yLimits));
stat=logical(Sig_front_ji);
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