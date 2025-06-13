function [corr_values, lags, max_corr, optimal_lag, trial_corrs] = electrode_diff_correlation(e1_dev, e1_std, e2_dev, e2_std, max_lag, srate, plot_results)
% ELECTRODE_DIFF_CORRELATION Calculates cross-correlation between the
% difference signals (deviant-standard) for a pair of electrodes using trial-level analysis
%
% INPUTS:
%   e1_dev       - Deviant signal from electrode 1 [trials x samples]
%   e1_std       - Standard signal from electrode 1 [trials x samples]
%   e2_dev       - Deviant signal from electrode 2 [trials x samples]
%   e2_std       - Standard signal from electrode 2 [trials x samples]
%   max_lag      - Maximum lag to compute in samples (default: 20)
%   srate        - Sampling rate in Hz (default: 1000)
%   plot_results - Whether to plot results (default: false)
%
% OUTPUTS:
%   corr_values - Cross-correlation values (averaged across trials)
%   lags        - Lag values in samples
%   max_corr    - Maximum correlation value
%   optimal_lag - Lag at maximum correlation (samples)
%   trial_corrs - Individual trial correlations [trials x lags]

% Handle optional arguments
if nargin < 5 || isempty(max_lag)
    max_lag = 20; % Default max lag (20 samples)
end

if nargin < 6 || isempty(srate)
    srate = 1000; % Default sampling rate (1000 Hz)
end

if nargin < 7 || isempty(plot_results)
    plot_results = false; % Default to no plotting for cluster jobs
end

% Check input dimensions
if size(e1_dev, 1) ~= size(e2_dev, 1) || size(e1_std, 1) ~= size(e2_std, 1)
    error('Number of trials must match between inputs');
end

% Get number of trials for each condition
n_dev_trials = size(e1_dev, 1);
n_std_trials = size(e1_std, 1);

% Determine how many trial pairs to use
if n_dev_trials == n_std_trials
    % Equal number of trials - use direct pairing
    n_pairs = n_dev_trials;
    pair_method = 'direct';
else
    % Unequal trials - use all possible combinations
    n_pairs = n_dev_trials * n_std_trials;
    pair_method = 'all_combinations';
end

% Initialize arrays
trial_corrs = [];

if strcmp(pair_method, 'direct')
    % Process each matched trial pair
    trial_corrs = zeros(n_pairs, 2*max_lag+1);
    
    for i = 1:n_pairs
        % Calculate difference signals for this trial pair
        diff_signal1 = e1_dev(i,:) - e1_std(i,:);
        diff_signal2 = e2_dev(i,:) - e2_std(i,:);
        
        % Compute cross-correlation
        [c, l] = xcorr(diff_signal1, diff_signal2, max_lag, 'coeff');
        trial_corrs(i,:) = c;
        
        % Store lags from first trial (same for all)
        if i == 1
            lags = l;
        end
    end
else
    % Process all possible combinations of deviant and standard trials
    trial_idx = 1;
    trial_corrs = zeros(n_pairs, 2*max_lag+1);
    
    for i = 1:n_dev_trials
        for j = 1:n_std_trials
            % Calculate difference signals
            diff_signal1 = e1_dev(i,:) - e1_std(j,:);
            diff_signal2 = e2_dev(i,:) - e2_std(j,:);
            
            % Compute cross-correlation
            [c, l] = xcorr(diff_signal1, diff_signal2, max_lag, 'coeff');
            trial_corrs(trial_idx,:) = c;
            
            % Store lags from first trial (same for all)
            if trial_idx == 1
                lags = l;
            end
            
            trial_idx = trial_idx + 1;
        end
    end
end

% Average correlations across all trial pairs
corr_values = mean(trial_corrs, 1);

% Find maximum correlation and its lag
[max_val, max_idx] = max(abs(corr_values));
max_corr = corr_values(max_idx); % Keep the actual sign
optimal_lag = lags(max_idx);

% Convert lag to time
lag_ms = optimal_lag * (1000/srate);

% Display results
fprintf('Maximum correlation: %.4f at lag of %d samples (%.1f ms)\n', ...
    max_corr, optimal_lag, lag_ms);
fprintf('Analysis performed on %d trial pairs\n', size(trial_corrs, 1));

% Skip plotting for cluster jobs unless explicitly requested
if plot_results
    figure('Position', [100 100 900 600]);
    
    % Main correlation plot
    subplot(2,1,1);
    hold on;
    
    % Calculate standard error of the mean
    sem = std(trial_corrs, 0, 1) / sqrt(size(trial_corrs, 1));
    
    % Add SEM shading
    x_region = [lags, fliplr(lags)];
    y_region_sem = [corr_values + sem, fliplr(corr_values - sem)];
    fill(x_region, y_region_sem, [0.7 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    % Plot the mean correlation
    plot(lags, corr_values, 'b-', 'LineWidth', 2);
    
    % Mark maximum correlation
    plot(optimal_lag, max_corr, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    % Add reference lines
    plot([0 0], ylim, 'k--');
    plot(xlim, [0 0], 'k--');
    
    % Add labels
    title(sprintf('Cross-correlation (max r=%.2f at lag=%d samples, %.1f ms)', ...
        max_corr, optimal_lag, lag_ms));
    xlabel('Lag (samples)');
    ylabel('Correlation');
    
    % Add legend
    legend({'SEM', 'Mean Correlation', 'Maximum'}, 'Location', 'Best');
    grid on;
    
    % Plot heatmap of individual trial correlations
    subplot(2,1,2);
    imagesc(lags, 1:size(trial_corrs,1), trial_corrs);
    colormap(redblue);
    c = colorbar;
    ylabel(c, 'Correlation');
    
    % Add reference line at zero lag
    hold on;
    plot([0 0], ylim, 'k--', 'LineWidth', 1);
    
    % Add labels
    title('Individual Trial Pair Correlations');
    xlabel('Lag (samples)');
    ylabel('Trial Pair');
    
    % Mark optimal lag
    plot([optimal_lag optimal_lag], ylim, 'r--', 'LineWidth', 1);
    
    % Add overall title
    sgtitle('Electrode Difference Signal Cross-Correlation Analysis');
end
end

% Helper function for colormap
function cmap = redblue(m)
% Creates a red-blue diverging colormap centered at zero
if nargin < 1
    m = 256;
end

% Define colors
red = [1 0 0];     % Red for positive values
blue = [0 0 1];    % Blue for negative values
white = [1 1 1];   % White for zero

% Create half maps
n = floor(m/2);

% Create gradient from blue to white
neg_map = zeros(n, 3);
for i = 1:3
    neg_map(:, i) = linspace(blue(i), white(i), n);
end

% Create gradient from white to red
pos_map = zeros(n, 3);
for i = 1:3
    pos_map(:, i) = linspace(white(i), red(i), n);
end

% Combine the maps
if mod(m, 2) == 1
    cmap = [neg_map; white; pos_map];
else
    cmap = [neg_map; pos_map];
end
end