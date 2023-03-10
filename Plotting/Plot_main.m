%% PLOTS COI/MI WITH PERMUTATION TESTING FOR SPECIFIED PARTICIPANT AND CONDITION
%...in the way we want to have it in the papers 

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


