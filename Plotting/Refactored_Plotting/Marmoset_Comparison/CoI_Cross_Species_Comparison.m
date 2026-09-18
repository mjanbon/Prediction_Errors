%% CoI Cross Series Comparison
% 
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(fileparts(this_dir)));
addpath(this_dir);
addpath(fileparts(this_dir));
addpath(fullfile(project_root, 'Code'));

%% SETUP AND PARAMETERS
% Load CoI comparison data with significance masks from permutation testing
USING_HPC = 0;  % 0 = local machine, 1 = Cambridge HPC, 2 = QMUL HPC
get_elec = 0;   % 0 = don't get electrodes of interest (not needed for plotting)
times = -25:74; % Time vector in ms (baseline: -25 to 0, post-stim: 0 to 74)
stim_idx = 26;  % Index of stimulus onset (time = 0 ms)
n_flies = 8;    % Total number of flies in the dataset
cutoff = 8;     % Cutoff for peripheral/central division

%% LOAD PARAMETERS
% Load fly parameters
[basefold, datatype, all_con, condition, subject, participants, EoI, ...
    srate, activity_tag, deviant_group_number, standard_group_number, corrected, stim_onset, baseline, ...
    start_cut_off, end_cut_off, kperm] = Max_get_param(USING_HPC, get_elec);
% Load marmoset parameters
csv_map = struct('peripheral', 'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\41467_2024_48329_MOESM4_ESM_Figure_3a_peripheral.csv',...
'central', 'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\41467_2024_48329_MOESM4_ESM_Figure_3c_central.csv');
% Output directory
out_dir = "C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\CoI-pipeline\Data\Results\Species_Comparison";

% Upscaling factor applied to Drosophila CoI matrices.
scale_factor = 4;

% Stim onset index for marmoset CoI matrix; only post-onset is used.
marmoset_stim_onset_ms = 100;

% If a marmoset CSV has been downsampled relative to the original
% millisecond/index grid, set the factor per region. For example, a 450 x
% 450 central matrix saved as 150 x 150 has central = 3, while the
% peripheral matrix can stay at 1 if it was not downsampled.
marmoset_downsample_factor = struct('peripheral', 1, 'central', 3);

% Resampling settings for permutation p-values and bootstrap confidence intervals.
Nperm = 50;
Nboot = 50;

%% RUN ANALYSIS
 [Cross_species_comparison_results] = compare_drosophila_marmoset_coI(participants,...
 basefold, datatype, activity_tag, condition, cutoff, times, csv_map, out_dir, scale_factor, marmoset_stim_onset_ms, marmoset_downsample_factor, Nperm, Nboot);

 %% PLOT MARMOSET COI
% plot_marmoset_central_coi_like_flies('C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\41467_2024_48329_MOESM4_ESM_Figure_3c_central.csv',basefold);
% plot_marmoset_central_coi_like_flies('C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\41467_2024_48329_MOESM4_ESM_Figure_3a_peripheral.csv',basefold,'Marmoset','Peripheral');
