%% set up color
function plot_line_CoI(basefold, datatype, participantnum, partname, condition, cutoff, times)
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

% These text and line sizes are best suited to subplot(4,4,x),
% which makes them about the right size for printing in a paper when saved to an A4 pdf.
[basefold, datatype, participantnum, ~, condition, participants, ~, re_epoch, dev_epochs, std_epochs, epoch_length, srate, low_cutoff, high_cutoff, filt_order, ~, start_cut_off, end_cut_off, kperm] = Get_param(0);

conditions = {'XY' 'XY_BB' };
for conditioni = 1:length(conditions)
for participanti = 1:length(participants)
[Highest] = Get_highest_MI(basefold,datatype,condition, participanti, char(participants(participanti)), cutoff, times);
HighAF.(char(participants(participanti))) = Highest;

[tempFFiall, tempFFi, temp_all_mask, tempFFm, tempFFmr, tempFFms,tempFFmi1, tempFFmi2,tempFFir,tempFFis, temp_nb, front_nb, tempfront_nb, tempint_nb, frontint_nb, ...
          frontFFiall, frontFFi, front_all_mask, frontFFm, frontFFmr, frontFFms, frontFFmi1, frontFFmi2, frontFFir, frontFFis,...
          tempfrontFFiall, tempfrontFFi, tempfrontFFm, tempfrontFFmr, tempfrontFFms, tempfrontFFmi1, tempfrontFFmi2, tempfrontFFir, tempfrontFFis,...
          tempintFFiall, tempintFFi, tempintFFm, tempintFFmr, tempintFFms, tempintFFmi1, tempintFFmi2, tempintFFir, tempintFFis,...
          frontintFFiall, frontintFFi, frontintFFm, frontintFFmr, frontintFFms, frontintFFmi1, frontintFFmi2, frontintFFir, frontintFFis] = get_plotting_CoI(basefold,datatype, participanti, char(participants(participanti)), char(conditions(conditioni)), cutoff, times);

Temporo_front.(char(participants(participanti))).(char(conditions(conditioni))).data = tempfrontFFi;
Temporo_front.(char(participants(participanti))).(char(conditions(conditioni))).sigMask = tempfrontFFm;
Temporo_front.(char(participants(participanti))).(char(conditions(conditioni))).sigMaskR = tempfrontFFmr;
Temporo_front.(char(participants(participanti))).(char(conditions(conditioni))).sigMaskS = tempfrontFFms;
Temporo_front.(char(participants(participanti))).(char(conditions(conditioni))).MI1 = tempFFmi1;
Temporo_front.(char(participants(participanti))).(char(conditions(conditioni))).MI2 = frontFFmi1;
Temporo_front.(char(participants(participanti))).(char(conditions(conditioni))).nbchan = tempfront_nb;

participants(participanti)
conditions(conditioni)
end
end

%ERPS
ERPJi = Temporo_front.Ji.XY.data;
ERPNr = Temporo_front.Nr.XY.data;

BBJi = Temporo_front.Ji.XY_BB.data
BBNr = Temporo_front.Nr.XY_BB.data;


times = 700:2:1148;

%% Average ERP
% H = figure;
figure('visible', 'off');
set(gcf,'renderer','Painters','Position', [10 10 500 1000]);
clim = [-0.005 0.005];
subplot(4,2,1)
allERP = cat(2,diag(ERPJi),diag(ERPNr));
stdshade_acj(squeeze(allERP)',0.2,'k',times);
xlim([700 1148]);  ylim([-0.025 0.025]);
set(gca,'ytick',[-0.025 0.025], 'yticklabel', [-0.025 0 0.025], 'xtick',[800 950 1100], 'xticklabel', [0 150 300]);

subplot(4,2,2)
allBB = cat(2,diag(BBJi),diag(BBNr));
stdshade_acj(squeeze(allBB)',0.2,'k',times);
xlim([700 1148]);  ylim([-0.025 0.025]);
set(gca,'ytick',[-0.025 0 0.025], 'yticklabel', [-0.025 0 0.025], 'xtick',[800 950 1100], 'xticklabel', [0 150 300]);

subplot(4,2,3)
a = diag(ERPJi); b = diag(ERPNr);
a(a<0)=0; b(b<0)=0;
allERP = cat(2,a,b);
stdshade_acj(squeeze(allERP)',0.2,'r', times);
xlim([700 1148]);  ylim([-0.025 0.025]);
set(gca,'ytick',[-0.025 0 0.025], 'yticklabel', [-0.025 0 0.025], 'xtick',[800 950 1100], 'xticklabel', [0 150 300]);

subplot(4,2,4)
a = diag(BBJi); b = diag(BBNr);
a(a<0)=0; b(b<0)=0;
allBB = cat(2,a,b);
stdshade_acj(squeeze(allBB)',0.2,'r', times);
xlim([700 1148]);  ylim([-0.025 0.025]);
set(gca,'ytick',[-0.025 0 0.025], 'yticklabel', [-0.025 0 0.025], 'xtick',[800 950 1100], 'xticklabel', [0 150 300]);

subplot(4,2,5)
a = diag(ERPJi); b = diag(ERPNr);
a(a>0)=0; b(b>0)=0;
allERP = cat(2,a,b);
stdshade_acj(squeeze(allERP)',0.2,'b', times);
xlim([700 1148]);  ylim([-0.025 0.025]);
set(gca,'ytick',[-0.025 0 0.025], 'yticklabel', [-0.025 0 0.025], 'xtick',[800 950 1100], 'xticklabel', [0 150 300]);

subplot(4,2,6)
a = diag(BBJi); b = diag(BBNr);
a(a>0)=0; b(b>0)=0;
allBB = cat(2,a,b);
stdshade_acj(squeeze(allBB)',0.2,'b', times);
xlim([700 1148]);  ylim([-0.025 0.025]);
set(gca,'ytick',[-0.025 0 0.025], 'yticklabel', [-0.025 0 0.025], 'xtick',[800 950 1100], 'xticklabel', [0 150 300]);

%% save
print(gcf, '-dpdf', '-bestfit', ['coI_summary_spatial_ERP_BB_TemporoFrontal_AVERAGE.pdf']);
end
