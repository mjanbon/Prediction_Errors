%% ERP Power Spectra and Coherence Analysis
% This script loads ERP data for a pair of electrodes, calculates power spectra
% and coherence, and displays results by frequency band using a data-driven approach.

%% Setup parameters
clear all; close all;

% Get parameters from Max_get_param
[basefold, datatype, all_con, condition, subject, participants, EoI, ...
 srate, activity_tag, deviant_group_number, standard_group_number, corrected, ...
 stim_onset, baseline, start_cut_off, end_cut_off, kperm] = Max_get_param(0, 0);

% Override specific parameters if needed
fly_name = 'R060721';  % Specify fly to analyze
electrode_pair = [1, 3];  % Analyze electrodes E1 and E3

% Analysis parameters for coherence calculation
window_length = 64;   % Reduced from 128 to 64
overlap_percent = 50;  % Overlap percentage

% Frequency analysis approach - use fine-grained frequency bins
freq_resolution = 1;   % Resolution in Hz
freq_max = 80;         % Maximum frequency to analyze
freq_edges = 0:freq_resolution:freq_max;
freq_centers = freq_edges(1:end-1) + freq_resolution/2;

%% Load data
fprintf('Loading data for %s, condition: %s, activity tag: %s\n', fly_name, condition, activity_tag);
data_file = fullfile(basefold, datatype, [fly_name, '_', condition, '.mat']);

% Check if file exists
if ~exist(data_file, 'file')
    error('Data file not found: %s', data_file);
end

% Load deviant and standard trials
[dvt, std] = load_trials_from_group_hyper(data_file, deviant_group_number, standard_group_number, corrected, srate);

fprintf('Loaded %d deviant trials and %d standard trials\n', dvt.trials, std.trials);
fprintf('Number of channels: %d\n', dvt.nbchan);

%% Extract data for selected electrodes
if max(electrode_pair) > dvt.nbchan
    error('Electrode pair [%d, %d] exceeds available channels (%d)', electrode_pair(1), electrode_pair(2), dvt.nbchan);
end

% Extract data for the two electrodes
e1_dev = squeeze(dvt.data(electrode_pair(1), :, :));  % Time x Trials
e2_dev = squeeze(dvt.data(electrode_pair(2), :, :));
e1_std = squeeze(std.data(electrode_pair(1), :, :));
e2_std = squeeze(std.data(electrode_pair(2), :, :));

% Transpose to get Trials x Time
e1_dev = e1_dev';
e2_dev = e2_dev';
e1_std = e1_std';
e2_std = e2_std';

% Get time vector
time = dvt.times;

% Print data dimensions for debugging
fprintf('Signal dimensions:\n');
fprintf('e1_dev: %d trials x %d samples\n', size(e1_dev, 1), size(e1_dev, 2));
fprintf('e2_dev: %d trials x %d samples\n', size(e2_dev, 1), size(e2_dev, 2));
fprintf('e1_std: %d trials x %d samples\n', size(e1_std, 1), size(e1_std, 2));
fprintf('e2_std: %d trials x %d samples\n', size(e2_std, 1), size(e2_std, 2));
fprintf('Window length: %d samples\n', window_length);

%% Calculate power spectra and coherence
fprintf('Calculating power spectra and coherence...\n');

% Calculate power spectra and coherence
[dev_pxx1, f] = calculate_power_spectrum(e1_dev, srate, window_length, overlap_percent);
[dev_pxx2, ~] = calculate_power_spectrum(e2_dev, srate, window_length, overlap_percent);
[std_pxx1, ~] = calculate_power_spectrum(e1_std, srate, window_length, overlap_percent);
[std_pxx2, ~] = calculate_power_spectrum(e2_std, srate, window_length, overlap_percent);

% Calculate coherence
[dev_coh, ~] = calculate_electrode_coherence(e1_dev, e2_dev, 1/srate, window_length, overlap_percent);
[std_coh, ~] = calculate_electrode_coherence(e1_std, e2_std, 1/srate, window_length, overlap_percent);

%% Data-driven frequency band analysis
% Create frequency bins for analysis
n_bins = length(freq_centers);
dev_coh_binned = zeros(n_bins, 1);
std_coh_binned = zeros(n_bins, 1);

% Calculate mean coherence in each bin
for i = 1:n_bins
    bin_min = freq_edges(i);
    bin_max = freq_edges(i+1);
    bin_idx = (f >= bin_min) & (f < bin_max);
    
    if any(bin_idx)
        dev_coh_binned(i) = mean(dev_coh(bin_idx));
        std_coh_binned(i) = mean(std_coh(bin_idx));
    end
end

% Find peaks in coherence spectrum (data-driven bands)
[peaks, locs] = findpeaks(dev_coh_binned, 'MinPeakHeight', 0.2, 'MinPeakDistance', 3);
peak_freqs = freq_centers(locs);

fprintf('Detected coherence peaks at frequencies: ');
fprintf('%.1f Hz, ', peak_freqs);
fprintf('\n');

% Also identify traditional bands for comparison
traditional_bands = struct('delta', [0.5 4], 'theta', [4 8], 'alpha', [8 13], ...
                         'beta', [13 30], 'gamma', [30 80]);
band_names = fieldnames(traditional_bands);

%% Visualization
% Plot 1: ERPs
figure('Position', [50 50 1200 800]);
subplot(2, 2, 1);
plot(time, mean(e1_dev, 1), 'b-', 'LineWidth', 1.5);
hold on;
plot(time, mean(e1_std, 1), 'r-', 'LineWidth', 1.5);
title(sprintf('Electrode E%d - ERP', electrode_pair(1)));
xlabel('Time (ms)');
ylabel('Amplitude (μV)');
legend('Deviant', 'Standard');
grid on;

subplot(2, 2, 2);
plot(time, mean(e2_dev, 1), 'b-', 'LineWidth', 1.5);
hold on;
plot(time, mean(e2_std, 1), 'r-', 'LineWidth', 1.5);
title(sprintf('Electrode E%d - ERP', electrode_pair(2)));
xlabel('Time (ms)');
ylabel('Amplitude (μV)');
legend('Deviant', 'Standard');
grid on;

% Plot 2: Power Spectra
subplot(2, 2, 3);
semilogy(f, dev_pxx1, 'b-', 'LineWidth', 1.5);
hold on;
semilogy(f, std_pxx1, 'r-', 'LineWidth', 1.5);
title(sprintf('Electrode E%d - Power Spectrum', electrode_pair(1)));
xlabel('Frequency (Hz)');
ylabel('Power (μV²/Hz)');
legend('Deviant', 'Standard');
grid on;
xlim([0 50]);

subplot(2, 2, 4);
semilogy(f, dev_pxx2, 'b-', 'LineWidth', 1.5);
hold on;
semilogy(f, std_pxx2, 'r-', 'LineWidth', 1.5);
title(sprintf('Electrode E%d - Power Spectrum', electrode_pair(2)));
xlabel('Frequency (Hz)');
ylabel('Power (μV²/Hz)');
legend('Deviant', 'Standard');
grid on;
xlim([0 50]);

% Plot 3: Coherence spectrum with detected peaks
figure('Position', [100 100 1200 600]);
subplot(1, 2, 1);
plot(f, dev_coh, 'b-', 'LineWidth', 1.5);
hold on;
plot(f, std_coh, 'r-', 'LineWidth', 1.5);

% Highlight traditional frequency bands
colors = lines(length(band_names));
for b = 1:length(band_names)
    band = band_names{b};
    band_range = traditional_bands.(band);
    x_band = [band_range(1), band_range(2), band_range(2), band_range(1)];
    y_band = [0, 0, 1, 1];
    patch(x_band, y_band, colors(b,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

% Mark detected peaks
if ~isempty(peaks)
    plot(peak_freqs, peaks, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    % Add text labels for peaks
    for i = 1:length(peaks)
        text(peak_freqs(i), peaks(i)+0.05, sprintf('%.1f Hz', peak_freqs(i)), ...
             'HorizontalAlignment', 'center');
    end
end

title(sprintf('Coherence between E%d and E%d', electrode_pair(1), electrode_pair(2)));
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend('Deviant', 'Standard');
grid on;
xlim([0 50]);
ylim([0 1]);

% Plot 4: Binned coherence spectrum
subplot(1, 2, 2);
bar(freq_centers, dev_coh_binned, 'FaceAlpha', 0.7);
hold on;
bar(freq_centers, std_coh_binned, 'FaceAlpha', 0.7);

% Add traditional band markers
for b = 1:length(band_names)
    band = band_names{b};
    band_range = traditional_bands.(band);
    x_mid = mean(band_range);
    plot([x_mid, x_mid], [0, 1], '--', 'Color', colors(b,:), 'LineWidth', 1.5);
    text(x_mid, 0.95, band, 'HorizontalAlignment', 'center', 'Color', colors(b,:), 'FontWeight', 'bold');
end

title('Binned Coherence Spectrum');
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend('Deviant', 'Standard');
grid on;
xlim([0 50]);
ylim([0 1]);

%% Calculate and display coherence in sliding frequency windows
window_size = 5; % Hz
step_size = 2.5; % Hz
window_starts = 0:step_size:(freq_max-window_size);
window_centers = window_starts + window_size/2;
n_windows = length(window_starts);

% Calculate mean coherence in sliding windows
dev_coh_windows = zeros(n_windows, 1);
std_coh_windows = zeros(n_windows, 1);

for i = 1:n_windows
    win_start = window_starts(i);
    win_end = win_start + window_size;
    win_idx = (f >= win_start) & (f < win_end);
    
    if any(win_idx)
        dev_coh_windows(i) = mean(dev_coh(win_idx));
        std_coh_windows(i) = mean(std_coh(win_idx));
    end
end

% Plot sliding window analysis
figure('Position', [100 700 1200 600]);

% Plot sliding window coherence
subplot(1, 2, 1);
plot(window_centers, dev_coh_windows, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(window_centers, std_coh_windows, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);

% Add traditional band markers
for b = 1:length(band_names)
    band = band_names{b};
    band_range = traditional_bands.(band);
    x_mid = mean(band_range);
    plot([x_mid, x_mid], [0, 1], '--', 'Color', colors(b,:), 'LineWidth', 1.5);
    text(x_mid, 0.95, band, 'HorizontalAlignment', 'center', 'Color', colors(b,:), 'FontWeight', 'bold');
end

title(sprintf('Sliding Window Coherence Analysis (%d Hz windows)', window_size));
xlabel('Window Center Frequency (Hz)');
ylabel('Mean Coherence');
legend('Deviant', 'Standard');
grid on;
xlim([0 50]);
ylim([0 1]);

% Plot coherence difference (deviant - standard)
subplot(1, 2, 2);
bar(window_centers, dev_coh_windows - std_coh_windows);

% Add band markers
for b = 1:length(band_names)
    band = band_names{b};
    band_range = traditional_bands.(band);
    x_mid = mean(band_range);
    plot([x_mid, x_mid], ylim, '--', 'Color', colors(b,:), 'LineWidth', 1.5);
    text(x_mid, 0.95*max(ylim), band, 'HorizontalAlignment', 'center', 'Color', colors(b,:), 'FontWeight', 'bold');
end

title('Coherence Difference (Deviant - Standard)');
xlabel('Window Center Frequency (Hz)');
ylabel('Difference in Coherence');
grid on;
xlim([0 50]);

%% Helper Functions

function [pxx, f] = calculate_power_spectrum(signal, fs, window_length, overlap_percent)
% CALCULATE_POWER_SPECTRUM
% Calculates power spectrum for a signal
%
% INPUTS:
%   signal        - Time series [trials x samples]
%   fs            - Sampling frequency in Hz
%   window_length - Length of window for spectral estimation in samples
%   overlap_percent - Percentage of overlap between windows
%
% OUTPUTS:
%   pxx           - Power spectrum
%   f             - Frequency vector in Hz

    % Handle input arguments
    if nargin < 3
        window_length = min(256, size(signal, 2));
    end
    if nargin < 4
        overlap_percent = 50;
    end
    
    % Ensure window length doesn't exceed signal length
    signal_length = size(signal, 2);
    if window_length > signal_length
        fprintf('Warning: Window length (%d) exceeds signal length (%d). Reducing window length.\n', ...
                window_length, signal_length);
        window_length = floor(signal_length/2) * 2; % Ensure even length
    end
    
    % Convert overlap from percentage to samples
    overlap = round(window_length * overlap_percent/100);
    
    % Create Hanning window
    win = hann(window_length);
    
    % Calculate NFFT
    nfft = max(256, 2^nextpow2(window_length));
    
    % Initialize power spectrum
    pxx = zeros(nfft/2+1, 1);
    
    % Process each trial
    num_trials = size(signal, 1);
    for trial = 1:num_trials
        % Get current trial data
        x = signal(trial, :);
        
        % Calculate periodogram using Welch's method
        [pxx_trial, f] = pwelch(x, win', overlap, nfft, fs);
        
        % Add to total
        pxx = pxx + pxx_trial;
    end
    
    % Average across trials
    pxx = pxx / num_trials;
end

function [coh, f] = calculate_electrode_coherence(signal1, signal2, dt, window_length, overlap_percent)
% CALCULATE_ELECTRODE_COHERENCE
% Calculates coherence between two electrode signals, following the same 
% approach as the field_field_coherence function in the Python code
%
% INPUTS:
%   signal1       - Time series from first electrode [trials x samples]
%   signal2       - Time series from second electrode [trials x samples]
%   dt            - Time step (1/sampling_rate) in seconds
%   window_length - Length of window for spectral estimation in samples
%   overlap_percent - Percentage of overlap between windows
%
% OUTPUTS:
%   coh           - Coherence spectrum (values between 0 and 1)
%   f             - Frequency vector in Hz

    % Handle input arguments
    if nargin < 4
        window_length = min(256, size(signal1, 2));
    end
    if nargin < 5
        overlap_percent = 50;
    end
    
    % Ensure window length doesn't exceed signal length
    signal_length = size(signal1, 2);
    if window_length > signal_length
        fprintf('Warning: Window length (%d) exceeds signal length (%d). Reducing window length.\n', ...
                window_length, signal_length);
        window_length = floor(signal_length/2) * 2; % Ensure even length
    end
    
    % Convert overlap from percentage to samples
    overlap = round(window_length * overlap_percent/100);
    
    % Get dimensions
    [T, N] = size(signal1);  % T = number of trials, N = number of samples
    
    % Initialize arrays for spectral estimates
    nfft = max(256, 2^nextpow2(window_length));
    SYY = zeros(nfft/2+1, 1);  % Auto-spectrum of y
    SXX = zeros(nfft/2+1, 1);  % Auto-spectrum of x
    SYX = zeros(nfft/2+1, 1);  % Cross-spectrum (complex)
    
    % Create Hanning window
    win = hann(window_length);
    
    % Loop over trials
    for k = 1:T
        % Get trial data
        x_trial = signal1(k, :);
        y_trial = signal2(k, :);
        
        % Process the signal in windows with overlap
        num_windows = floor((length(x_trial) - overlap) / (window_length - overlap));
        
        % For each window
        for w = 1:num_windows
            % Extract window of data
            start_idx = (w-1) * (window_length - overlap) + 1;
            end_idx = start_idx + window_length - 1;
            
            if end_idx > length(x_trial)
                break;
            end
            
            % Extract window
            x_win = x_trial(start_idx:end_idx);
            y_win = y_trial(start_idx:end_idx);
            
            % Remove mean
            x_win = x_win - mean(x_win);
            y_win = y_win - mean(y_win);
            
            % Apply window
            x_win = x_win .* win';
            y_win = y_win .* win';
            
            % Calculate FFT
            X = fft(x_win, nfft);
            Y = fft(y_win, nfft);
            
            % Only use first half of spectrum (up to Nyquist frequency)
            X = X(1:nfft/2+1);
            Y = Y(1:nfft/2+1);
            
            % Accumulate cross-spectrum and auto-spectra
            SYX = SYX + (Y .* conj(X));
            SXX = SXX + (X .* conj(X));
            SYY = SYY + (Y .* conj(Y));
        end
    end
    
    % Normalize by number of trials and windows
    total_segments = T * num_windows;
    SYX = SYX / total_segments;
    SXX = SXX / total_segments;
    SYY = SYY / total_segments;
    
    % Calculate coherence
    coh = abs(SYX).^2 ./ (SXX .* SYY);
    
    % Create frequency vector
    f = (0:nfft/2) / (N * dt);
end