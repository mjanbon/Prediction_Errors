function analyze_synergy_redundancy_across_flies(fly_list, basefold, datatype, activity_tag, condition, timing)
% ANALYZE_SYNERGY_REDUNDANCY_ACROSS_FLIES
% Analyzes how synergy and redundancy in CoI matrices change with electrode distance
% across multiple flies, and plots results with error bars.
%
% INPUTS:
%   fly_list         - Cell array of fly names to include
%   basefold         - Base folder path 
%   datatype         - Type of data
%   activity_tag     - Activity tag
%   condition        - Experimental condition
%   timing           - Time vector

    fprintf('Analyzing synergy and redundancy by electrode distance across %d flies\n', length(fly_list));
    
    % Initialize aggregation structures
    max_distance = 20; % Maximum electrode distance to analyze
    distances = 1:max_distance;
    
    % For each fly, we'll track these metrics by distance
    synergy_by_fly = zeros(length(fly_list), max_distance);
    redundancy_by_fly = zeros(length(fly_list), max_distance);
    peak_synergy_by_fly = zeros(length(fly_list), max_distance);
    peak_redundancy_by_fly = zeros(length(fly_list), max_distance);
    synergy_percent_by_fly = zeros(length(fly_list), max_distance);
    redundancy_percent_by_fly = zeros(length(fly_list), max_distance);
    count_pairs_by_fly = zeros(length(fly_list), max_distance);
    % Also collect per-pair (per-electrode) metrics pooled across flies
    all_total_syn = cell(1, max_distance);
    all_total_red = cell(1, max_distance);
    all_peak_syn = cell(1, max_distance);
    all_peak_red = cell(1, max_distance);
    
    % Process each fly
    valid_fly_count = 0;
    
    for f = 1:length(fly_list)
        participant_name = fly_list{f};
        try
            fprintf('Processing fly %s (%d of %d)\n', participant_name, f, length(fly_list));
            
            % Load CoI data
            [~, CoI] = max_load_data(basefold, datatype, condition, activity_tag, participant_name);
            
            % Load electrodes of interest (EoI)
            eoi_file = fullfile(basefold, 'DataEOI', ['EoI_data_', datatype, '.mat']);
            load(eoi_file, 'EoI');
            
            if ~isfield(EoI, participant_name) || ~isfield(EoI.(participant_name), activity_tag) || ...
               ~isfield(EoI.(participant_name).(activity_tag), condition)
                fprintf('  Skipping fly %s - no electrode data found\n', participant_name);
                continue;
            end
            
            all_elecs = EoI.(participant_name).(activity_tag).(condition);
            
            % Extract electrode numbers and sort them
            elec_nums = zeros(length(all_elecs), 1);
            for i = 1:length(all_elecs)
                elec_name = char(all_elecs(i));
                % Extract number from electrode name (assuming format like 'E1', 'E2', etc.)
                matches = regexp(elec_name, '\d+', 'match');
                if ~isempty(matches)
                    elec_nums(i) = str2double(matches{1});
                else
                    elec_nums(i) = i; % Default to index if no number found
                end
            end
            
            % Sort electrodes by their numbers
            [sorted_nums, sort_idx] = sort(elec_nums);
            sorted_elecs = all_elecs(sort_idx);
            
            % Initialize data structures for this fly
            n_elecs = length(sorted_elecs);
            fly_max_distance = min(n_elecs - 1, max_distance);
            
            % Initialize matrix storage and metrics for this fly
            avg_matrix_by_distance = cell(fly_max_distance, 1);
            count_by_distance = zeros(fly_max_distance, 1);
            
            % Extract full CoI matrix for each electrode pair by distance
            for i = 1:n_elecs
                for j = i+1:n_elecs
                    % Calculate distance between electrodes
                    distance = abs(sorted_nums(j) - sorted_nums(i));
                    if distance > fly_max_distance
                        continue; % Skip distances beyond our analysis range
                    end
                    
                    E1 = char(sorted_elecs(i));
                    E2 = char(sorted_elecs(j));
                    
                    % Try both electrode name orderings
                    pair_names = {[E1, '_', E2], [E2, '_', E1]};
                    
                    for p = 1:length(pair_names)
                        elec_pair = pair_names{p};
                        if isfield(CoI.(participant_name).(condition).data, elec_pair)
                            % Get CoI matrix - full time-time matrix
                            coi_matrix = CoI.(participant_name).(condition).data.(elec_pair);
                            % Initialize the average matrix for this distance if first pair
                            if count_by_distance(distance) == 0
                                avg_matrix_by_distance{distance} = zeros(size(coi_matrix));
                            end

                            % Add to average matrix (per-fly average)
                            avg_matrix_by_distance{distance} = avg_matrix_by_distance{distance} + coi_matrix;
                            count_by_distance(distance) = count_by_distance(distance) + 1;

                            % ALSO record this pair's metrics in pooled lists (use raw matrix)
                            try
                                syn_pair_sum = sum(coi_matrix(coi_matrix < 0));
                                red_pair_sum = sum(coi_matrix(coi_matrix > 0));
                                peak_syn_pair = min(coi_matrix(:));
                                peak_red_pair = max(coi_matrix(:));
                            catch
                                syn_pair_sum = NaN; red_pair_sum = NaN; peak_syn_pair = NaN; peak_red_pair = NaN;
                            end
                            if isempty(all_total_syn{distance}), all_total_syn{distance} = syn_pair_sum; else all_total_syn{distance}(end+1) = syn_pair_sum; end
                            if isempty(all_total_red{distance}), all_total_red{distance} = red_pair_sum; else all_total_red{distance}(end+1) = red_pair_sum; end
                            if isempty(all_peak_syn{distance}), all_peak_syn{distance} = peak_syn_pair; else all_peak_syn{distance}(end+1) = peak_syn_pair; end
                            if isempty(all_peak_red{distance}), all_peak_red{distance} = peak_red_pair; else all_peak_red{distance}(end+1) = peak_red_pair; end
                            break; % Found the pair, so break the inner loop
                        end
                    end
                end
            end
            
            % Calculate metrics for each distance for this fly
            for d = 1:fly_max_distance
                if count_by_distance(d) > 0
                    % Calculate average matrix
                    avg_matrix_by_distance{d} = avg_matrix_by_distance{d} / count_by_distance(d);
                    
                    % Synergy matrix (only negative values)
                    synergy_matrix = avg_matrix_by_distance{d} .* (avg_matrix_by_distance{d} < 0);
                    
                    % Redundancy matrix (only positive values)
                    redundancy_matrix = avg_matrix_by_distance{d} .* (avg_matrix_by_distance{d} > 0);
                    
                    % Total synergy (sum of negative values)
                    synergy_by_fly(f, d) = sum(synergy_matrix(:));
                    
                    % Total redundancy (sum of positive values)
                    redundancy_by_fly(f, d) = sum(redundancy_matrix(:));
                    
                    % Peak synergy (most negative value)
                    [peak_synergy_by_fly(f, d), ~] = min(avg_matrix_by_distance{d}(:));
                    
                    % Peak redundancy (most positive value)
                    [peak_redundancy_by_fly(f, d), ~] = max(avg_matrix_by_distance{d}(:));
                    
                    % Percentage of matrix showing synergy
                    synergy_percent_by_fly(f, d) = 100 * sum(synergy_matrix(:) < 0) / numel(synergy_matrix);
                    
                    % Percentage of matrix showing redundancy
                    redundancy_percent_by_fly(f, d) = 100 * sum(redundancy_matrix(:) > 0) / numel(redundancy_matrix);
                    
                    % Store count of pairs at this distance
                    count_pairs_by_fly(f, d) = count_by_distance(d);
                end
            end
            
            valid_fly_count = valid_fly_count + 1;
            
        catch e
            fprintf('Error processing fly %s: %s\n', participant_name, e.message);
            disp(getReport(e, 'extended'));
        end
    end
    
    if valid_fly_count == 0
        error('No valid flies processed');
    end
    
    % Calculate statistics across flies
    synergy_mean = mean(synergy_by_fly, 1, 'omitnan');
    % SEM across flies (keeps previous behaviour)
    synergy_sem = std(synergy_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(synergy_by_fly), 1));
    % Pooled SEM using all electrode pairs (for shaded error regions)
    pooled_total_syn_sem = nan(1, max_distance);
    pooled_total_red_sem = nan(1, max_distance);
    pooled_peak_syn_sem = nan(1, max_distance);
    pooled_peak_red_sem = nan(1, max_distance);
    for d = 1:max_distance
        vals = all_total_syn{d};
        if ~isempty(vals)
            pooled_total_syn_sem(d) = std(vals, 0) / sqrt(numel(vals));
        else
            pooled_total_syn_sem(d) = NaN;
        end
        vals = all_total_red{d};
        if ~isempty(vals)
            pooled_total_red_sem(d) = std(vals, 0) / sqrt(numel(vals));
        else
            pooled_total_red_sem(d) = NaN;
        end
        vals = all_peak_syn{d};
        if ~isempty(vals)
            pooled_peak_syn_sem(d) = std(vals, 0) / sqrt(numel(vals));
        else
            pooled_peak_syn_sem(d) = NaN;
        end
        vals = all_peak_red{d};
        if ~isempty(vals)
            pooled_peak_red_sem(d) = std(vals, 0) / sqrt(numel(vals));
        else
            pooled_peak_red_sem(d) = NaN;
        end
    end
    
    redundancy_mean = mean(redundancy_by_fly, 1, 'omitnan');
    redundancy_sem = std(redundancy_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(redundancy_by_fly), 1));
    
    peak_synergy_mean = mean(peak_synergy_by_fly, 1, 'omitnan');
    peak_synergy_sem = std(peak_synergy_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(peak_synergy_by_fly), 1));
    
    peak_redundancy_mean = mean(peak_redundancy_by_fly, 1, 'omitnan');
    peak_redundancy_sem = std(peak_redundancy_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(peak_redundancy_by_fly), 1));
    
    syn_percent_mean = mean(synergy_percent_by_fly, 1, 'omitnan');
    syn_percent_sem = std(synergy_percent_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(synergy_percent_by_fly), 1));
    
    red_percent_mean = mean(redundancy_percent_by_fly, 1, 'omitnan');
    red_percent_sem = std(redundancy_percent_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(redundancy_percent_by_fly), 1));
    
    count_mean = mean(count_pairs_by_fly, 1, 'omitnan');
    
    % Create and save figures
    plot_synergy_redundancy_results(distances, synergy_mean, synergy_sem, redundancy_mean, redundancy_sem, ...
                                  peak_synergy_mean, peak_synergy_sem, peak_redundancy_mean, peak_redundancy_sem,...
                                  syn_percent_mean, syn_percent_sem, red_percent_mean, red_percent_sem, ... 
                                  count_mean, pooled_total_syn_sem, pooled_total_red_sem, pooled_peak_syn_sem, pooled_peak_red_sem, ...
                                  activity_tag, condition, basefold, valid_fly_count);

    % Report and save exponential decay fits (mean curves) and exponential mixed model
    save_dir = fullfile(basefold, 'Results', 'CoI_distance', 'Across_electrodes');
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end

    try
        % Exponential fit on means: y = a * exp(b * distance)
        dvec = distances(:);
        valid_d = count_mean(:) > 0;
        dfit = dvec(valid_d);
        y_syn_mean = -synergy_mean(valid_d)';
        y_red_mean = redundancy_mean(valid_d)';

        exp_fit_summary = struct();

        valid_syn = isfinite(y_syn_mean) & (y_syn_mean > 0);
        if sum(valid_syn) >= 2
            mdl_exp_syn = fitlm(dfit(valid_syn), log(y_syn_mean(valid_syn)));
            ci_syn = coefCI(mdl_exp_syn);
            b_syn = mdl_exp_syn.Coefficients.Estimate(2);
            b_syn_ci = ci_syn(2,:);
            exp_fit_summary.syn = struct( ...
                'a', exp(mdl_exp_syn.Coefficients.Estimate(1)), ...
                'b', b_syn, ...
                'b_ci', b_syn_ci, ...
                'p', mdl_exp_syn.Coefficients.pValue(2), ...
                'R2', mdl_exp_syn.Rsquared.Ordinary, ...
                'fold_per_distance', exp(b_syn), ...
                'fold_ci', exp(b_syn_ci));
            fprintf('Exponential decay fit (synergy): b=%g, 95%% CI [%g, %g], p=%g, R2=%g\n', ...
                b_syn, b_syn_ci(1), b_syn_ci(2), mdl_exp_syn.Coefficients.pValue(2), mdl_exp_syn.Rsquared.Ordinary);
        else
            mdl_exp_syn = [];
            exp_fit_summary.syn = struct('a', NaN, 'b', NaN, 'b_ci', [NaN NaN], 'p', NaN, 'R2', NaN, 'fold_per_distance', NaN, 'fold_ci', [NaN NaN]);
            fprintf('Exponential decay fit (synergy): insufficient valid points\n');
        end

        valid_red = isfinite(y_red_mean) & (y_red_mean > 0);
        if sum(valid_red) >= 2
            mdl_exp_red = fitlm(dfit(valid_red), log(y_red_mean(valid_red)));
            ci_red = coefCI(mdl_exp_red);
            b_red = mdl_exp_red.Coefficients.Estimate(2);
            b_red_ci = ci_red(2,:);
            exp_fit_summary.red = struct( ...
                'a', exp(mdl_exp_red.Coefficients.Estimate(1)), ...
                'b', b_red, ...
                'b_ci', b_red_ci, ...
                'p', mdl_exp_red.Coefficients.pValue(2), ...
                'R2', mdl_exp_red.Rsquared.Ordinary, ...
                'fold_per_distance', exp(b_red), ...
                'fold_ci', exp(b_red_ci));
            fprintf('Exponential decay fit (redundancy): b=%g, 95%% CI [%g, %g], p=%g, R2=%g\n', ...
                b_red, b_red_ci(1), b_red_ci(2), mdl_exp_red.Coefficients.pValue(2), mdl_exp_red.Rsquared.Ordinary);
        else
            mdl_exp_red = [];
            exp_fit_summary.red = struct('a', NaN, 'b', NaN, 'b_ci', [NaN NaN], 'p', NaN, 'R2', NaN, 'fold_per_distance', NaN, 'fold_ci', [NaN NaN]);
            fprintf('Exponential decay fit (redundancy): insufficient valid points\n');
        end

        save(fullfile(save_dir, sprintf('%s_%s_ExponentialDecayFits_%dFlies.mat', activity_tag, condition, valid_fly_count)), ...
            'exp_fit_summary', 'mdl_exp_syn', 'mdl_exp_red');
    catch ME
        warning('analyze_synergy_redundancy:expFitReportFailed', 'Failed to report/save exponential fits: %s', ME.message);
    end

    try
        % Exponential mixed model via log-linear LME on fly-level totals
        Y_syn = -synergy_by_fly;   % make positive
        Y_red = redundancy_by_fly;
        Fly = {};
        Distance = [];
        MetricType = {};
        Value = [];
        for f = 1:length(fly_list)
            for d = 1:max_distance
                if isfinite(Y_syn(f,d)) && (Y_syn(f,d) > 0)
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Distance(end+1,1) = d; %#ok<SAGROW>
                    MetricType{end+1,1} = 'syn'; %#ok<SAGROW>
                    Value(end+1,1) = Y_syn(f,d); %#ok<SAGROW>
                end
                if isfinite(Y_red(f,d)) && (Y_red(f,d) > 0)
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Distance(end+1,1) = d; %#ok<SAGROW>
                    MetricType{end+1,1} = 'red'; %#ok<SAGROW>
                    Value(end+1,1) = Y_red(f,d); %#ok<SAGROW>
                end
            end
        end

        Texp = table(categorical(Fly), Distance(:), categorical(MetricType), Value(:), ...
            'VariableNames', {'Fly','Distance','MetricType','Value'});
        Texp.MetricType = categorical(Texp.MetricType, {'syn','red'});
        Texp.LogValue = log(Texp.Value);

        % Interaction tests whether log-decay slope differs between metrics
        lme_exp = fitlme(Texp, 'LogValue ~ Distance*MetricType + (1|Fly)');
        disp('Exponential mixed model fixed effects (log scale):');
        disp(lme_exp.Coefficients);

        coeffsTbl = lme_exp.Coefficients;
        names = coeffsTbl.Name;
        idxD = find(strcmp(names, 'Distance'));
        idxInt = find(strcmp(names, 'Distance:MetricType_red'));
        betaVec = lme_exp.Coefficients.Estimate;
        covB = lme_exp.CoefficientCovariance;

        if ~isempty(idxD)
            slope_syn_log = betaVec(idxD);
            var_syn_log = covB(idxD, idxD);
            se_syn_log = sqrt(var_syn_log);
            ci_syn_log = slope_syn_log + [-1, 1] * 1.96 * se_syn_log;
        else
            slope_syn_log = NaN; se_syn_log = NaN; ci_syn_log = [NaN NaN];
        end

        if ~isempty(idxInt) && ~isempty(idxD)
            slope_red_log = slope_syn_log + betaVec(idxInt);
            var_red_log = covB(idxD,idxD) + covB(idxInt,idxInt) + 2 * covB(idxD,idxInt);
            se_red_log = sqrt(var_red_log);
            ci_red_log = slope_red_log + [-1, 1] * 1.96 * se_red_log;
            interaction_est = betaVec(idxInt);
            interaction_se = coeffsTbl.SE(idxInt);
            interaction_p = coeffsTbl.pValue(idxInt);
            interaction_ci = interaction_est + [-1, 1] * 1.96 * interaction_se;
        else
            slope_red_log = NaN; se_red_log = NaN; ci_red_log = [NaN NaN];
            interaction_est = NaN; interaction_se = NaN; interaction_p = NaN; interaction_ci = [NaN NaN];
        end

        exp_mixed_summary = struct();
        exp_mixed_summary.slope_syn_log = slope_syn_log;
        exp_mixed_summary.se_syn_log = se_syn_log;
        exp_mixed_summary.ci_syn_log = ci_syn_log;
        exp_mixed_summary.fold_syn = exp(slope_syn_log);
        exp_mixed_summary.fold_syn_ci = exp(ci_syn_log);
        exp_mixed_summary.slope_red_log = slope_red_log;
        exp_mixed_summary.se_red_log = se_red_log;
        exp_mixed_summary.ci_red_log = ci_red_log;
        exp_mixed_summary.fold_red = exp(slope_red_log);
        exp_mixed_summary.fold_red_ci = exp(ci_red_log);
        exp_mixed_summary.interaction_est = interaction_est;
        exp_mixed_summary.interaction_se = interaction_se;
        exp_mixed_summary.interaction_ci = interaction_ci;
        exp_mixed_summary.interaction_p = interaction_p;

        fprintf('Exponential mixed model slopes (log scale): syn=%g (95%% CI [%g, %g]); red=%g (95%% CI [%g, %g])\n', ...
            slope_syn_log, ci_syn_log(1), ci_syn_log(2), slope_red_log, ci_red_log(1), ci_red_log(2));
        fprintf('Distance x MetricType(red) interaction: est=%g, SE=%g, 95%% CI [%g, %g], p=%g\n', ...
            interaction_est, interaction_se, interaction_ci(1), interaction_ci(2), interaction_p);

        save(fullfile(save_dir, sprintf('%s_%s_ExponentialMixedModel_%dFlies.mat', activity_tag, condition, valid_fly_count)), ...
            'lme_exp', 'Texp', 'exp_mixed_summary');
    catch ME
        warning('analyze_synergy_redundancy:expMixedFailed', 'Exponential mixed model failed: %s', ME.message);
    end

    try
        % Exponential fit on PEAK means: y = a * exp(b * distance)
        dvec = distances(:);
        valid_d = count_mean(:) > 0;
        dfit = dvec(valid_d);
        y_peak_syn_mean = -peak_synergy_mean(valid_d)';
        y_peak_red_mean = peak_redundancy_mean(valid_d)';

        exp_fit_peak_summary = struct();

        valid_peak_syn = isfinite(y_peak_syn_mean) & (y_peak_syn_mean > 0);
        if sum(valid_peak_syn) >= 2
            mdl_exp_peak_syn = fitlm(dfit(valid_peak_syn), log(y_peak_syn_mean(valid_peak_syn)));
            ci_peak_syn = coefCI(mdl_exp_peak_syn);
            b_peak_syn = mdl_exp_peak_syn.Coefficients.Estimate(2);
            b_peak_syn_ci = ci_peak_syn(2,:);
            exp_fit_peak_summary.syn = struct( ...
                'a', exp(mdl_exp_peak_syn.Coefficients.Estimate(1)), ...
                'b', b_peak_syn, ...
                'b_ci', b_peak_syn_ci, ...
                'p', mdl_exp_peak_syn.Coefficients.pValue(2), ...
                'R2', mdl_exp_peak_syn.Rsquared.Ordinary, ...
                'fold_per_distance', exp(b_peak_syn), ...
                'fold_ci', exp(b_peak_syn_ci));
            fprintf('Exponential decay fit (peak synergy): b=%g, 95%% CI [%g, %g], p=%g, R2=%g\n', ...
                b_peak_syn, b_peak_syn_ci(1), b_peak_syn_ci(2), mdl_exp_peak_syn.Coefficients.pValue(2), mdl_exp_peak_syn.Rsquared.Ordinary);
        else
            mdl_exp_peak_syn = [];
            exp_fit_peak_summary.syn = struct('a', NaN, 'b', NaN, 'b_ci', [NaN NaN], 'p', NaN, 'R2', NaN, 'fold_per_distance', NaN, 'fold_ci', [NaN NaN]);
            fprintf('Exponential decay fit (peak synergy): insufficient valid points\n');
        end

        valid_peak_red = isfinite(y_peak_red_mean) & (y_peak_red_mean > 0);
        if sum(valid_peak_red) >= 2
            mdl_exp_peak_red = fitlm(dfit(valid_peak_red), log(y_peak_red_mean(valid_peak_red)));
            ci_peak_red = coefCI(mdl_exp_peak_red);
            b_peak_red = mdl_exp_peak_red.Coefficients.Estimate(2);
            b_peak_red_ci = ci_peak_red(2,:);
            exp_fit_peak_summary.red = struct( ...
                'a', exp(mdl_exp_peak_red.Coefficients.Estimate(1)), ...
                'b', b_peak_red, ...
                'b_ci', b_peak_red_ci, ...
                'p', mdl_exp_peak_red.Coefficients.pValue(2), ...
                'R2', mdl_exp_peak_red.Rsquared.Ordinary, ...
                'fold_per_distance', exp(b_peak_red), ...
                'fold_ci', exp(b_peak_red_ci));
            fprintf('Exponential decay fit (peak redundancy): b=%g, 95%% CI [%g, %g], p=%g, R2=%g\n', ...
                b_peak_red, b_peak_red_ci(1), b_peak_red_ci(2), mdl_exp_peak_red.Coefficients.pValue(2), mdl_exp_peak_red.Rsquared.Ordinary);
        else
            mdl_exp_peak_red = [];
            exp_fit_peak_summary.red = struct('a', NaN, 'b', NaN, 'b_ci', [NaN NaN], 'p', NaN, 'R2', NaN, 'fold_per_distance', NaN, 'fold_ci', [NaN NaN]);
            fprintf('Exponential decay fit (peak redundancy): insufficient valid points\n');
        end

        save(fullfile(save_dir, sprintf('%s_%s_PeakExponentialDecayFits_%dFlies.mat', activity_tag, condition, valid_fly_count)), ...
            'exp_fit_peak_summary', 'mdl_exp_peak_syn', 'mdl_exp_peak_red');
    catch ME
        warning('analyze_synergy_redundancy:expPeakFitReportFailed', 'Failed to report/save peak exponential fits: %s', ME.message);
    end

    try
        % Exponential mixed model for PEAK values via log-linear LME
        Y_peak_syn = -peak_synergy_by_fly;   % make positive
        Y_peak_red = peak_redundancy_by_fly;
        Fly = {};
        Distance = [];
        MetricType = {};
        Value = [];
        for f = 1:length(fly_list)
            for d = 1:max_distance
                if isfinite(Y_peak_syn(f,d)) && (Y_peak_syn(f,d) > 0)
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Distance(end+1,1) = d; %#ok<SAGROW>
                    MetricType{end+1,1} = 'syn'; %#ok<SAGROW>
                    Value(end+1,1) = Y_peak_syn(f,d); %#ok<SAGROW>
                end
                if isfinite(Y_peak_red(f,d)) && (Y_peak_red(f,d) > 0)
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Distance(end+1,1) = d; %#ok<SAGROW>
                    MetricType{end+1,1} = 'red'; %#ok<SAGROW>
                    Value(end+1,1) = Y_peak_red(f,d); %#ok<SAGROW>
                end
            end
        end

        Texp_peak = table(categorical(Fly), Distance(:), categorical(MetricType), Value(:), ...
            'VariableNames', {'Fly','Distance','MetricType','Value'});
        Texp_peak.MetricType = categorical(Texp_peak.MetricType, {'syn','red'});
        Texp_peak.LogValue = log(Texp_peak.Value);

        lme_exp_peak = fitlme(Texp_peak, 'LogValue ~ Distance*MetricType + (1|Fly)');
        disp('Exponential mixed model (peak) fixed effects (log scale):');
        disp(lme_exp_peak.Coefficients);

        coeffsTbl = lme_exp_peak.Coefficients;
        names = coeffsTbl.Name;
        idxD = find(strcmp(names, 'Distance'));
        idxInt = find(strcmp(names, 'Distance:MetricType_red'));
        betaVec = lme_exp_peak.Coefficients.Estimate;
        covB = lme_exp_peak.CoefficientCovariance;

        if ~isempty(idxD)
            slope_syn_log_peak = betaVec(idxD);
            var_syn_log_peak = covB(idxD, idxD);
            se_syn_log_peak = sqrt(var_syn_log_peak);
            ci_syn_log_peak = slope_syn_log_peak + [-1, 1] * 1.96 * se_syn_log_peak;
        else
            slope_syn_log_peak = NaN; se_syn_log_peak = NaN; ci_syn_log_peak = [NaN NaN];
        end

        if ~isempty(idxInt) && ~isempty(idxD)
            slope_red_log_peak = slope_syn_log_peak + betaVec(idxInt);
            var_red_log_peak = covB(idxD,idxD) + covB(idxInt,idxInt) + 2 * covB(idxD,idxInt);
            se_red_log_peak = sqrt(var_red_log_peak);
            ci_red_log_peak = slope_red_log_peak + [-1, 1] * 1.96 * se_red_log_peak;
            interaction_est_peak = betaVec(idxInt);
            interaction_se_peak = coeffsTbl.SE(idxInt);
            interaction_p_peak = coeffsTbl.pValue(idxInt);
            interaction_ci_peak = interaction_est_peak + [-1, 1] * 1.96 * interaction_se_peak;
        else
            slope_red_log_peak = NaN; se_red_log_peak = NaN; ci_red_log_peak = [NaN NaN];
            interaction_est_peak = NaN; interaction_se_peak = NaN; interaction_p_peak = NaN; interaction_ci_peak = [NaN NaN];
        end

        exp_mixed_peak_summary = struct();
        exp_mixed_peak_summary.slope_syn_log = slope_syn_log_peak;
        exp_mixed_peak_summary.se_syn_log = se_syn_log_peak;
        exp_mixed_peak_summary.ci_syn_log = ci_syn_log_peak;
        exp_mixed_peak_summary.fold_syn = exp(slope_syn_log_peak);
        exp_mixed_peak_summary.fold_syn_ci = exp(ci_syn_log_peak);
        exp_mixed_peak_summary.slope_red_log = slope_red_log_peak;
        exp_mixed_peak_summary.se_red_log = se_red_log_peak;
        exp_mixed_peak_summary.ci_red_log = ci_red_log_peak;
        exp_mixed_peak_summary.fold_red = exp(slope_red_log_peak);
        exp_mixed_peak_summary.fold_red_ci = exp(ci_red_log_peak);
        exp_mixed_peak_summary.interaction_est = interaction_est_peak;
        exp_mixed_peak_summary.interaction_se = interaction_se_peak;
        exp_mixed_peak_summary.interaction_ci = interaction_ci_peak;
        exp_mixed_peak_summary.interaction_p = interaction_p_peak;

        fprintf('Exponential mixed model (peak) slopes (log scale): syn=%g (95%% CI [%g, %g]); red=%g (95%% CI [%g, %g])\n', ...
            slope_syn_log_peak, ci_syn_log_peak(1), ci_syn_log_peak(2), slope_red_log_peak, ci_red_log_peak(1), ci_red_log_peak(2));
        fprintf('Peak distance x MetricType(red) interaction: est=%g, SE=%g, 95%% CI [%g, %g], p=%g\n', ...
            interaction_est_peak, interaction_se_peak, interaction_ci_peak(1), interaction_ci_peak(2), interaction_p_peak);

        save(fullfile(save_dir, sprintf('%s_%s_PeakExponentialMixedModel_%dFlies.mat', activity_tag, condition, valid_fly_count)), ...
            'lme_exp_peak', 'Texp_peak', 'exp_mixed_peak_summary');
    catch ME
        warning('analyze_synergy_redundancy:expPeakMixedFailed', 'Peak exponential mixed model failed: %s', ME.message);
    end
    
    %% Optional: Violin plots for redundancy by region (only if redundancy_all provided)
    if exist('redundancy_all','var') && iscell(redundancy_all)
        [n_regions_rows, ~] = size(redundancy_all);
        n_regions = n_regions_rows;

        % Default region labels and activity tags if not provided
        region_labels = arrayfun(@(i)sprintf('Region %d', i), 1:n_regions, 'UniformOutput', false)';
        activity_tags = {'Wake', 'Sleep'};

        % Ensure save directory exists (reuse same Results folder)
        save_dir = fullfile(basefold, 'Results', 'CoI_distance', 'Across_electrodes');
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end

        for r = 1:n_regions
            region_label = region_labels{r};
            data_wake = redundancy_all{r,1};
            data_sleep = redundancy_all{r,2};

            % Prepare cell array for violin
            data_cell = {data_wake, data_sleep};

            figure;
            violin(data_cell, 'xlabel', activity_tags);
            ylabel('Redundancy (FFi > 0)');
            title(['Redundancy for ', region_label, ' (all pairs, all flies)']);

            % Unpaired t-test
            try
                [~,p] = ttest2(data_wake, data_sleep);
            catch
                p = NaN;
            end
            fprintf('%s: unpaired t-test p=%.4g\n', region_label, p);

            % Annotate plot
            y_max = max([data_wake(:); data_sleep(:)]);
            if isempty(y_max) || ~isfinite(y_max)
                y_max = 0;
            end
            text(1.5, y_max, sprintf('p=%.3g', p), 'HorizontalAlignment', 'center');

            % Save figure
            saveas(gcf, fullfile(save_dir, ['Redundancy_Violin_', region_label, '.png']));
            saveas(gcf, fullfile(save_dir, ['Redundancy_Violin_', region_label, '.fig']));
            print(gcf, fullfile(save_dir, ['Redundancy_Violin_', region_label, '.svg']), '-dsvg');
            close(gcf);
        end
    else
        % No redundancy_all provided; skip violin plots gracefully
    end
end

function plot_synergy_redundancy_results(distances, synergy_mean, synergy_sem, redundancy_mean, redundancy_sem, ... 
                                       peak_synergy_mean, peak_synergy_sem, peak_redundancy_mean, peak_redundancy_sem, ...
                                       ~, ~, ~, ~, ... 
                                       count_mean, pooled_total_syn_sem, pooled_total_red_sem, pooled_peak_syn_sem, pooled_peak_red_sem, ...
                                       activity_tag, condition, basefold, n_flies)
    % Find the maximum distance that has data
    valid_distances = find(count_mean > 0);
    if isempty(valid_distances)
        error('No valid distance data found');
    end
    max_valid_distance = max(valid_distances);
    
    % Create x-axis limited to valid distances
    plot_distances = distances(1:max_valid_distance);
    
    % Plot 1: Synergy and Redundancy together
    fig1 = figure('Position', [100 100 1200 800]);
    
    % Plot 1: Total Synergy and Redundancy by distance
    subplot(2, 2, 1);
    hold on;
    
    % Synergy with shaded error region (negative values, so multiply by -1 for plotting)
    h_syn = plot(plot_distances, -synergy_mean(1:max_valid_distance), 'b-', 'LineWidth', 0.3);
    hold on;
    % Use pooled SEM (all electrode pairs) for shading when available
    shade_error_region(plot_distances, -synergy_mean(1:max_valid_distance), pooled_total_syn_sem(1:max_valid_distance), [0 0.4 0.8]);
    
    % Redundancy with shaded error region
    h_red = plot(plot_distances, redundancy_mean(1:max_valid_distance), 'r-', 'LineWidth', 0.3);
    shade_error_region(plot_distances, redundancy_mean(1:max_valid_distance), pooled_total_red_sem(1:max_valid_distance), [0.8 0.2 0.2]);
    
    % Fit linear and logarithmic models to the change over distance
    % We'll model absolute synergy (plot uses -synergy_mean) and redundancy
    try
        y_syn = -synergy_mean(1:max_valid_distance); % make positive for fitting/visualisation
        y_red = redundancy_mean(1:max_valid_distance);
        x = plot_distances(:);

        % Linear fits
        mdl_lin_syn = fitlm(x, y_syn(:));
        mdl_lin_red = fitlm(x, y_red(:));

        % Exponential decay fits: model y = a * exp(b*x).
        % Fit by linearizing: log(y) = log(a) + b*x. Only fit where y>0.
        valid_syn = isfinite(y_syn) & (y_syn > 0);
        valid_red = isfinite(y_red) & (y_red > 0);
        if sum(valid_syn) >= 2
            mdl_exp_syn = fitlm(x(valid_syn), log(y_syn(valid_syn)));
        else
            mdl_exp_syn = [];
        end
        if sum(valid_red) >= 2
            mdl_exp_red = fitlm(x(valid_red), log(y_red(valid_red)));
        else
            mdl_exp_red = [];
        end

        % Generate smooth x for plotting fitted curves
        xfit = linspace(min(x), max(x), 200)';
        y_lin_syn_fit = predict(mdl_lin_syn, xfit);
        if ~isempty(mdl_exp_syn)
            y_log_syn_fit = exp(predict(mdl_exp_syn, xfit));
        else
            y_log_syn_fit = nan(size(xfit));
        end
        y_lin_red_fit = predict(mdl_lin_red, xfit);
        if ~isempty(mdl_exp_red)
            y_log_red_fit = exp(predict(mdl_exp_red, xfit));
        else
            y_log_red_fit = nan(size(xfit));
        end

    % Do not draw linear/exponential fit lines on the figure.

        % Prepare fit statistics to save
        fit_stats = struct();
        fit_stats.lin_syn = struct('slope', mdl_lin_syn.Coefficients.Estimate(2), 'p', mdl_lin_syn.Coefficients.pValue(2), 'R2', mdl_lin_syn.Rsquared.Ordinary);
        if ~isempty(mdl_exp_syn)
            fit_stats.exp_syn = struct('a', exp(mdl_exp_syn.Coefficients.Estimate(1)), 'b', mdl_exp_syn.Coefficients.Estimate(2), 'p', mdl_exp_syn.Coefficients.pValue(2), 'R2', mdl_exp_syn.Rsquared.Ordinary);
        else
            fit_stats.exp_syn = struct('a', NaN, 'b', NaN, 'p', NaN, 'R2', NaN);
        end
        fit_stats.lin_red = struct('slope', mdl_lin_red.Coefficients.Estimate(2), 'p', mdl_lin_red.Coefficients.pValue(2), 'R2', mdl_lin_red.Rsquared.Ordinary);
        if ~isempty(mdl_exp_red)
            fit_stats.exp_red = struct('a', exp(mdl_exp_red.Coefficients.Estimate(1)), 'b', mdl_exp_red.Coefficients.Estimate(2), 'p', mdl_exp_red.Coefficients.pValue(2), 'R2', mdl_exp_red.Rsquared.Ordinary);
        else
            fit_stats.exp_red = struct('a', NaN, 'b', NaN, 'p', NaN, 'R2', NaN);
        end

        % (Removed) fit summary annotation box - omitted to keep figure clean

        % Save fit stats
        save_dir = fullfile(basefold, 'Results', 'CoI_distance', 'Across_electrodes');
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        save(fullfile(save_dir, sprintf('%s_%s_SynRed_fits_%dFlies.mat', activity_tag, condition, n_flies)), 'fit_stats');
    catch ME
        warning('analyze_synergy_redundancy:fitFailed', 'Could not compute/plot fits: %s', ME.message);
    end

    % Create a compact legend showing only the main lines (synergy, redundancy)
    % Place legend inside the axes (northeast) so it does not shift subplot positions
    legend([h_syn, h_red], {'Synergy','Redundancy'}, 'Location', 'northeast', 'FontSize', 9, 'Box', 'off');
    xlabel('Number of separating electrodes');
    ylabel('Total Information (bits)');
    title('Total Synergy & Redundancy by Distance');
    grid off;
    xlim([0.5, max_valid_distance+0.5]);
    
    % Plot 2: Peak Synergy and Redundancy by distance
    subplot(2, 2, 2);
    hold on;
    
    % Peak Synergy with shaded error region
    h_ps = plot(plot_distances, -peak_synergy_mean(1:max_valid_distance), 'b-', 'LineWidth', 0.3); %#ok<NASGU>
    shade_error_region(plot_distances, -peak_synergy_mean(1:max_valid_distance), pooled_peak_syn_sem(1:max_valid_distance), [0 0.4 0.8]);

    % Peak Redundancy with shaded error region
    h_pr = plot(plot_distances, peak_redundancy_mean(1:max_valid_distance), 'r-', 'LineWidth', 0.3); %#ok<NASGU>
    shade_error_region(plot_distances, peak_redundancy_mean(1:max_valid_distance), pooled_peak_red_sem(1:max_valid_distance), [0.8 0.2 0.2]);

    % Do not create legend for peak subplot in combined figure (avoids 'data*' labels)
    xlabel('Number of separating electrodes');
    ylabel('Peak Value (bits)');
    title('Peak Synergy & Redundancy by Distance');
    grid off;
    xlim([0.5, max_valid_distance+0.5]);

    % Also fit linear and exponential models to the peak values and overlay
    try
        y_peak_syn = -peak_synergy_mean(1:max_valid_distance);
        y_peak_red = peak_redundancy_mean(1:max_valid_distance);
        x = plot_distances(:);

        % Linear fits for peaks
        mdl_lin_peak_syn = fitlm(x, y_peak_syn(:));
        mdl_lin_peak_red = fitlm(x, y_peak_red(:));

        % Exponential fits for peaks (log-linear) where y>0
        valid_ps = isfinite(y_peak_syn) & (y_peak_syn > 0);
        valid_pr = isfinite(y_peak_red) & (y_peak_red > 0);
        if sum(valid_ps) >= 2
            mdl_exp_peak_syn = fitlm(x(valid_ps), log(y_peak_syn(valid_ps)));
        else
            mdl_exp_peak_syn = [];
        end
        if sum(valid_pr) >= 2
            mdl_exp_peak_red = fitlm(x(valid_pr), log(y_peak_red(valid_pr)));
        else
            mdl_exp_peak_red = [];
        end

        % Predict for smooth curves
        xfit = linspace(min(x), max(x), 200)';
        y_lin_peak_syn_fit = predict(mdl_lin_peak_syn, xfit);
        y_lin_peak_red_fit = predict(mdl_lin_peak_red, xfit);
        if ~isempty(mdl_exp_peak_syn)
            y_exp_peak_syn_fit = exp(predict(mdl_exp_peak_syn, xfit));
        else
            y_exp_peak_syn_fit = nan(size(xfit));
        end
        if ~isempty(mdl_exp_peak_red)
            y_exp_peak_red_fit = exp(predict(mdl_exp_peak_red, xfit));
        else
            y_exp_peak_red_fit = nan(size(xfit));
        end

        % Do not draw linear/exponential fit lines on the figure.

        % Add peak fit stats to fit_stats if it exists
        if exist('fit_stats','var')
            fit_stats.lin_peak_syn = struct('slope', mdl_lin_peak_syn.Coefficients.Estimate(2), 'p', mdl_lin_peak_syn.Coefficients.pValue(2), 'R2', mdl_lin_peak_syn.Rsquared.Ordinary);
            fit_stats.lin_peak_red = struct('slope', mdl_lin_peak_red.Coefficients.Estimate(2), 'p', mdl_lin_peak_red.Coefficients.pValue(2), 'R2', mdl_lin_peak_red.Rsquared.Ordinary);
            if ~isempty(mdl_exp_peak_syn)
                fit_stats.exp_peak_syn = struct('a', exp(mdl_exp_peak_syn.Coefficients.Estimate(1)), 'b', mdl_exp_peak_syn.Coefficients.Estimate(2), 'p', mdl_exp_peak_syn.Coefficients.pValue(2), 'R2', mdl_exp_peak_syn.Rsquared.Ordinary);
            else
                fit_stats.exp_peak_syn = struct('a', NaN, 'b', NaN, 'p', NaN, 'R2', NaN);
            end
            if ~isempty(mdl_exp_peak_red)
                fit_stats.exp_peak_red = struct('a', exp(mdl_exp_peak_red.Coefficients.Estimate(1)), 'b', mdl_exp_peak_red.Coefficients.Estimate(2), 'p', mdl_exp_peak_red.Coefficients.pValue(2), 'R2', mdl_exp_peak_red.Rsquared.Ordinary);
            else
                fit_stats.exp_peak_red = struct('a', NaN, 'b', NaN, 'p', NaN, 'R2', NaN);
            end
            % Save updated fit_stats if save_dir exists
            if exist('save_dir','var') && ~isempty(save_dir)
                save(fullfile(save_dir, sprintf('%s_%s_SynRed_fits_%dFlies.mat', activity_tag, condition, n_flies)), 'fit_stats');
            end
        end
    catch ME
        warning('analyze_synergy_redundancy:peakFitFailed', '%s', ME.message);
    end
    
    % (Removed) Percentage and pair-count subplots: omitted per requested change
    
    % Global title
    sgtitle(sprintf('%s %s: Synergy & Redundancy Analysis by Distance (n=%d flies)', activity_tag, condition, n_flies));
    
    % Plot 2: Separate plots for synergy metrics (percentage panel removed)
    fig2 = figure('Position', [100 300 1000 400]);

    % Plot 1: Total Synergy by distance
    subplot(1, 2, 1);
    hold on;
    h_tot_syn = plot(plot_distances, -synergy_mean(1:max_valid_distance), 'b-', 'LineWidth', 0.3); %#ok<NASGU>
    shade_error_region(plot_distances, -synergy_mean(1:max_valid_distance), synergy_sem(1:max_valid_distance), [0 0.4 0.8]);
    xlabel('Number of separating electrodes');
    ylabel('Total Synergy (bits)');
    title('Total Synergy by Distance');
    grid off;
    xlim([0.5, max_valid_distance+0.5]);

    % Plot 2: Peak Synergy by distance
    subplot(1, 2, 2);
    hold on;
    h_peak_syn = plot(plot_distances, -peak_synergy_mean(1:max_valid_distance), 'b-', 'LineWidth', 0.3); %#ok<NASGU>
    shade_error_region(plot_distances, -peak_synergy_mean(1:max_valid_distance), peak_synergy_sem(1:max_valid_distance), [0 0.4 0.8]);
    xlabel('Number of separating electrodes');
    ylabel('Peak Synergy (bits)');
    title('Peak Synergy by Distance');
    grid off;
    xlim([0.5, max_valid_distance+0.5]);

    % No legend for separate synergy figure (user requested no legends)

    % Global title for synergy (percentage panel removed)
    sgtitle(sprintf('%s %s: Synergy Analysis by Distance (n=%d flies)', activity_tag, condition, n_flies));
    
    % Plot 3: Separate plots for redundancy metrics (percentage panel removed)
    fig3 = figure('Position', [100 600 1000 400]);

    % Plot 1: Total Redundancy by distance
    subplot(1, 2, 1);
    hold on;
    h_tot_red = plot(plot_distances, redundancy_mean(1:max_valid_distance), 'r-', 'LineWidth', 0.3); %#ok<NASGU>
    shade_error_region(plot_distances, redundancy_mean(1:max_valid_distance), redundancy_sem(1:max_valid_distance), [0.8 0.2 0.2]);
    xlabel('Number of separating electrodes');
    ylabel('Total Redundancy (bits)');
    title('Total Redundancy by Distance');
    grid off;
    xlim([0.5, max_valid_distance+0.5]);

    % Plot 2: Peak Redundancy by distance
    subplot(1, 2, 2);
    hold on;
    h_peak_red = plot(plot_distances, peak_redundancy_mean(1:max_valid_distance), 'r-', 'LineWidth', 0.3); %#ok<NASGU>
    shade_error_region(plot_distances, peak_redundancy_mean(1:max_valid_distance), peak_redundancy_sem(1:max_valid_distance), [0.8 0.2 0.2]);
    xlabel('Number of separating electrodes');
    ylabel('Peak Redundancy (bits)');
    title('Peak Redundancy by Distance');
    grid off;
    xlim([0.5, max_valid_distance+0.5]);

    % No legend for separate redundancy figure (user requested no legends)

    % Global title for redundancy (percentage panel removed)
    sgtitle(sprintf('%s %s: Redundancy Analysis by Distance (n=%d flies)', activity_tag, condition, n_flies));
    
    % Save figures
    save_dir = fullfile(basefold, 'Results', 'CoI_distance', 'Across_electrodes');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    set(findall([fig1 fig2 fig3], 'Type', 'Line'), 'LineWidth', 0.3);
    set(findall([fig1 fig2 fig3], 'Type', 'Axes'), 'XGrid', 'off', 'YGrid', 'off', 'LineWidth', 0.25);

    % Save combined figure
    saveas(fig1, fullfile(save_dir, sprintf('%s_%s_SynRedByDistance_Combined_%dFlies.fig', activity_tag, condition, n_flies)), 'fig');
    saveas(fig1, fullfile(save_dir, sprintf('%s_%s_SynRedByDistance_Combined_%dFlies.png', activity_tag, condition, n_flies)), 'png');
    print(fig1, fullfile(save_dir, sprintf('%s_%s_SynRedByDistance_Combined_%dFlies.svg', activity_tag, condition, n_flies)), '-dsvg', '-painters');
    
    % Save synergy figure
    saveas(fig2, fullfile(save_dir, sprintf('%s_%s_SynergyByDistance_%dFlies.fig', activity_tag, condition, n_flies)), 'fig');
    saveas(fig2, fullfile(save_dir, sprintf('%s_%s_SynergyByDistance_%dFlies.png', activity_tag, condition, n_flies)), 'png');
    print(fig2, fullfile(save_dir, sprintf('%s_%s_SynergyByDistance_%dFlies.svg', activity_tag, condition, n_flies)), '-dsvg', '-painters');
    
    % Save redundancy figure
    saveas(fig3, fullfile(save_dir, sprintf('%s_%s_RedundancyByDistance_%dFlies.fig', activity_tag, condition, n_flies)), 'fig');
    saveas(fig3, fullfile(save_dir, sprintf('%s_%s_RedundancyByDistance_%dFlies.png', activity_tag, condition, n_flies)), 'png');
    print(fig3, fullfile(save_dir, sprintf('%s_%s_RedundancyByDistance_%dFlies.svg', activity_tag, condition, n_flies)), '-dsvg', '-painters');
    
    fprintf('Saved figures to %s\n', save_dir);
end

function shade_error_region(x, y, err, color)
    % Create x coordinates for the error region (forward and back)
    x_polygon = [x, fliplr(x)];
    
    % Create y coordinates for the error region (upper and lower bounds)
    y_upper = y + err;
    y_lower = y - err;
    y_polygon = [y_upper, fliplr(y_lower)];
    
    % Plot the error region as a filled polygon
    h = fill(x_polygon, y_polygon, color);
    set(h, 'EdgeColor', 'none');
    set(h, 'FaceAlpha', 0.3); % Transparency level
    % Prevent the filled polygon from appearing in legends
    try
        set(h, 'HandleVisibility', 'off');
        lh = get(h, 'Annotation');
        if isstruct(lh) || ~isempty(lh)
            set(get(h, 'Annotation').LegendInformation, 'IconDisplayStyle', 'off');
        end
    catch
        % Ignore if properties aren't available in older MATLAB versions
    end
end

function analyze_synergy_redundancy_by_electrode(fly_list, basefold, datatype, activity_tag, condition)
% ANALYZE_SYNERGY_REDUNDANCY_BY_ELECTRODE
% For each electrode (E1..En) compute the average CoI matrix across all
% pairs that include that electrode, then compute total and peak synergy
% (negative values) and redundancy (positive values). Aggregate across
% flies, plot the metrics across electrode positions (retina -> central: Emax->E1),
% save figures, CSV and fit statistics.

    fprintf('Analyzing synergy and redundancy by electrode across %d flies\n', length(fly_list));

    % First pass: discover maximum electrode number across flies to build a
    % consistent electrode axis (handles cases where some flies have fewer electrodes).
    max_elec_num = 0;
    EoI_file = fullfile(basefold, 'DataEOI', ['EoI_data_', datatype, '.mat']);
    if exist(EoI_file, 'file')
        load(EoI_file, 'EoI');
    else
        error('EoI file not found: %s', EoI_file);
    end

    for f = 1:length(fly_list)
        name = fly_list{f};
        if isfield(EoI, name) && isfield(EoI.(name), activity_tag) && isfield(EoI.(name).(activity_tag), condition)
            all_elecs = EoI.(name).(activity_tag).(condition);
            for i = 1:length(all_elecs)
                matches = regexp(char(all_elecs(i)), '\\d+', 'match');
                if ~isempty(matches)
                    max_elec_num = max(max_elec_num, str2double(matches{1}));
                end
            end
        end
    end

    if max_elec_num == 0
        error('No electrode numbers discovered across flies');
    end

    % Preallocate arrays: flies x electrode_number (1..max_elec_num)
    n_flies = length(fly_list);
    total_syn_by_fly = nan(n_flies, max_elec_num);
    total_red_by_fly = nan(n_flies, max_elec_num);
    peak_syn_by_fly = nan(n_flies, max_elec_num);
    peak_red_by_fly = nan(n_flies, max_elec_num);
    count_pairs_by_fly = zeros(n_flies, max_elec_num);

    % Iterate flies and accumulate per-electrode average matrices
    for f = 1:n_flies
        participant_name = fly_list{f};
        try
            fprintf('  Processing %s (%d/%d)\n', participant_name, f, n_flies);
            [~, CoI] = max_load_data(basefold, datatype, condition, activity_tag, participant_name);

            if ~isfield(EoI, participant_name) || ~isfield(EoI.(participant_name), activity_tag) || ...
               ~isfield(EoI.(participant_name).(activity_tag), condition)
                fprintf('    No EoI for %s - skipping\n', participant_name);
                continue;
            end

            elec_list = EoI.(participant_name).(activity_tag).(condition);
            % Map electrode names to numbers and keep their string names
            elec_nums = nan(length(elec_list),1);
            elec_names = cell(length(elec_list),1);
            for i = 1:length(elec_list)
                s = char(elec_list(i));
                elec_names{i} = s;
                m = regexp(s, '\\d+', 'match');
                if ~isempty(m)
                    elec_nums(i) = str2double(m{1});
                else
                    elec_nums(i) = i;
                end
            end

            % For each electrode number present for this fly, find all pairs
            % in the CoI data that include that electrode and average their
            % full time-time CoI matrices.
            for ei = 1:length(elec_list)
                ename = char(elec_names{ei});
                enum = elec_nums(ei);
                if isnan(enum)
                    continue;
                end

                accum = [];
                cnt = 0;

                % Search for self-pair field in CoI data: only use matrices where
                % the electrode is paired with itself (e.g., 'E1_E1'). This
                % produces one CoI matrix per electrode per fly.
                fields = fieldnames(CoI.(participant_name).(condition).data);
                for pf = 1:length(fields)
                    fname = fields{pf};
                    parts = strsplit(fname, '_');
                    % Only accept exact self-pair fields like 'E3_E3'
                    if numel(parts) == 2 && strcmp(parts{1}, ename) && strcmp(parts{2}, ename)
                        mat = CoI.(participant_name).(condition).data.(fname);
                        if isempty(accum)
                            accum = zeros(size(mat));
                        end
                        accum = accum + mat;
                        cnt = cnt + 1;
                    end
                end

                if cnt > 0
                    avg_mat = accum ./ cnt;
                    syn_mat = avg_mat .* (avg_mat < 0);
                    red_mat = avg_mat .* (avg_mat > 0);

                    total_syn_by_fly(f, enum) = sum(syn_mat(:));
                    total_red_by_fly(f, enum) = sum(red_mat(:));
                    peak_syn_by_fly(f, enum) = min(avg_mat(:));
                    peak_red_by_fly(f, enum) = max(avg_mat(:));
                    count_pairs_by_fly(f, enum) = cnt;
                else
                    % leave as NaN / zero count
                end
            end

        catch E
            fprintf('    Error for %s: %s\n', participant_name, E.message);
            disp(getReport(E,'basic'));
        end
    end

    % Aggregate across flies (omit NaNs)
    mean_total_syn = nanmean(total_syn_by_fly, 1);
    sem_total_syn = nanstd(total_syn_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(total_syn_by_fly),1));
    mean_total_red = nanmean(total_red_by_fly, 1);
    sem_total_red = nanstd(total_red_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(total_red_by_fly),1));

    mean_peak_syn = nanmean(peak_syn_by_fly, 1);
    sem_peak_syn = nanstd(peak_syn_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(peak_syn_by_fly),1));
    mean_peak_red = nanmean(peak_red_by_fly, 1);
    sem_peak_red = nanstd(peak_red_by_fly, 1, 'omitnan') ./ sqrt(sum(~isnan(peak_red_by_fly),1));

    % Prepare x-axis: electrode numbers 1..max_elec_num but we want to plot
    % retina->central: Emax -> E1, so reverse for display
    elec_x = 1:max_elec_num;
    elec_x_plot = fliplr(elec_x);

    % Plot totals (combined)
    save_dir = fullfile(basefold, 'Results', 'CoI_distance', 'By_electrode');
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end

    fig = figure('Position', [100 100 1100 600]);
    subplot(1,2,1);
    hold on;
    % plot total synergy (make positive for plotting)
    y_syn_plot = -mean_total_syn(elec_x_plot);
    y_syn_sem = fliplr(sem_total_syn(elec_x_plot));
    y_red_plot = mean_total_red(elec_x_plot);
    y_red_sem = fliplr(sem_total_red(elec_x_plot));

    h1 = plot(elec_x, y_syn_plot, 'b-', 'LineWidth', 0.3); %#ok<NASGU>
    shade_error_region(elec_x, y_syn_plot, y_syn_sem, [0 0.4 0.8]);
    h2 = plot(elec_x, y_red_plot, 'r-', 'LineWidth', 0.3); %#ok<NASGU>
    shade_error_region(elec_x, y_red_plot, y_red_sem, [0.8 0.2 0.2]);

    xlabel('Electrode number (E_{retina} -> E_{central})');
    xticks(elec_x);
    xticklabels(arrayfun(@(n)sprintf('E%d', n), fliplr(elec_x), 'UniformOutput', false));
    xlim([min(elec_x)-0.5, max(elec_x)+0.5]);
    ylabel('Total Information (bits)');
    title(sprintf('%s %s: Total Synergy & Redundancy by Electrode (n=%d flies)', activity_tag, condition, n_flies));
    legend({'Synergy','Redundancy'}, 'Location', 'northeast', 'Box', 'off');
    grid off;

    % Fit linear and exponential trends across electrode number for totals
    try
        x_fit = elec_x';
        y_syn_fit = y_syn_plot(:);
        y_red_fit = y_red_plot(:);
        valid_syn = isfinite(y_syn_fit) & (y_syn_fit > 0);
        valid_red = isfinite(y_red_fit) & (y_red_fit > 0);

        mdl_lin_syn = fitlm(x_fit(valid_syn), y_syn_fit(valid_syn));
        mdl_lin_red = fitlm(x_fit(valid_red), y_red_fit(valid_red));
        x_smooth = linspace(min(x_fit), max(x_fit), 200)';
        plot(x_smooth, predict(mdl_lin_syn, x_smooth), '--', 'Color', [0 0.4 0.8]);
        plot(x_smooth, predict(mdl_lin_red, x_smooth), '--', 'Color', [0.8 0.2 0.2]);

        if sum(valid_syn) >= 2
            mdl_exp_syn = fitlm(x_fit(valid_syn), log(y_syn_fit(valid_syn)));
            plot(x_smooth, exp(predict(mdl_exp_syn, x_smooth)), ':', 'Color', [0 0.4 0.8]);
        else
            mdl_exp_syn = [];
        end
        if sum(valid_red) >=2
            mdl_exp_red = fitlm(x_fit(valid_red), log(y_red_fit(valid_red)));
            plot(x_smooth, exp(predict(mdl_exp_red, x_smooth)), ':', 'Color', [0.8 0.2 0.2]);
        else
            mdl_exp_red = [];
        end

        fit_stats = struct();
        fit_stats.lin_total_syn = struct('slope', mdl_lin_syn.Coefficients.Estimate(2), 'p', mdl_lin_syn.Coefficients.pValue(2), 'R2', mdl_lin_syn.Rsquared.Ordinary);
        fit_stats.lin_total_red = struct('slope', mdl_lin_red.Coefficients.Estimate(2), 'p', mdl_lin_red.Coefficients.pValue(2), 'R2', mdl_lin_red.Rsquared.Ordinary);
        if ~isempty(mdl_exp_syn)
            fit_stats.exp_total_syn = struct('a', exp(mdl_exp_syn.Coefficients.Estimate(1)), 'b', mdl_exp_syn.Coefficients.Estimate(2), 'p', mdl_exp_syn.Coefficients.pValue(2), 'R2', mdl_exp_syn.Rsquared.Ordinary);
        else
            fit_stats.exp_total_syn = struct('a', NaN, 'b', NaN, 'p', NaN, 'R2', NaN);
        end
        if ~isempty(mdl_exp_red)
            fit_stats.exp_total_red = struct('a', exp(mdl_exp_red.Coefficients.Estimate(1)), 'b', mdl_exp_red.Coefficients.Estimate(2), 'p', mdl_exp_red.Coefficients.pValue(2), 'R2', mdl_exp_red.Rsquared.Ordinary);
        else
            fit_stats.exp_total_red = struct('a', NaN, 'b', NaN, 'p', NaN, 'R2', NaN);
        end
    catch E
        warning('analyze_by_elec:fitFailed', '%s', E.message);
        fit_stats = struct();
    end

    % Plot peaks in second panel
    subplot(1,2,2);
    hold on;
    y_ps_plot = -fliplr(mean_peak_syn(elec_x));
    y_ps_sem = fliplr(sem_peak_syn(elec_x));
    y_pr_plot = fliplr(mean_peak_red(elec_x));
    y_pr_sem = fliplr(sem_peak_red(elec_x));

    plot(elec_x, y_ps_plot, 'b-', 'LineWidth', 0.3);
    shade_error_region(elec_x, y_ps_plot, y_ps_sem, [0 0.4 0.8]);
    plot(elec_x, y_pr_plot, 'r-', 'LineWidth', 0.3);
    shade_error_region(elec_x, y_pr_plot, y_pr_sem, [0.8 0.2 0.2]);
    xlabel('Electrode number (E_{retina} -> E_{central})');
    xticks(elec_x);
    xticklabels(arrayfun(@(n)sprintf('E%d', n), fliplr(elec_x), 'UniformOutput', false));
    xlim([min(elec_x)-0.5, max(elec_x)+0.5]);
    ylabel('Peak Value (bits)');
    title('Peak Synergy & Redundancy by Electrode');
    grid off;

    sgtitle(sprintf('%s %s: Synergy/Redundancy by Electrode (n=%d flies)', activity_tag, condition, n_flies));

    % Save outputs: figure, fits and CSV with per-fly per-electrode metrics
    set(findall(fig, 'Type', 'Line'), 'LineWidth', 0.3);
    set(findall(fig, 'Type', 'Axes'), 'XGrid', 'off', 'YGrid', 'off', 'LineWidth', 0.25);
    saveas(fig, fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_Combined_%dFlies.fig', activity_tag, condition, n_flies)));
    saveas(fig, fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_Combined_%dFlies.png', activity_tag, condition, n_flies)));
    print(fig, fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_Combined_%dFlies.svg', activity_tag, condition, n_flies)), '-dsvg', '-painters');
    save(fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_fits_%dFlies.mat', activity_tag, condition, n_flies)), 'fit_stats');

    % Build CSV rows
    csv_rows = {};
    header = {'Fly','Electrode','ElectrodeNumber','TotalSynergy','TotalRedundancy','PeakSynergy','PeakRedundancy','PairCount'};
    csv_rows{1} = header;
    rowi = 2;
    for f = 1:n_flies
        for e = 1:max_elec_num
            if ~isnan(total_syn_by_fly(f,e)) || ~isnan(total_red_by_fly(f,e)) || count_pairs_by_fly(f,e)>0
                csv_rows{rowi,1} = fly_list{f};
                csv_rows{rowi,2} = sprintf('E%d', e);
                csv_rows{rowi,3} = e;
                csv_rows{rowi,4} = total_syn_by_fly(f,e);
                csv_rows{rowi,5} = total_red_by_fly(f,e);
                csv_rows{rowi,6} = peak_syn_by_fly(f,e);
                csv_rows{rowi,7} = peak_red_by_fly(f,e);
                csv_rows{rowi,8} = count_pairs_by_fly(f,e);
                rowi = rowi + 1;
            end
        end
    end

    % Write CSV using low-level file IO to avoid table dependencies
    csv_file = fullfile(save_dir, sprintf('%s_%s_SynRed_ByElectrode.csv', activity_tag, condition));
    fid = fopen(csv_file, 'w');
    for r = 1:size(csv_rows,1)
        if r == 1
            fprintf(fid, '%s,', csv_rows{r,1:end-1});
            fprintf(fid, '%s\n', csv_rows{r,end});
        else
            % some cells might be numeric
            row = csv_rows(r,:);
            for c = 1:length(row)
                val = row{c};
                if isnumeric(val)
                    if isnan(val)
                        fprintf(fid, ',');
                    else
                        fprintf(fid, '%g,', val);
                    end
                else
                    % string
                    fprintf(fid, '%s,', val);
                end
            end
            % replace last comma with newline
            fseek(fid, -1, 'cof');
            fprintf(fid, '\n');
        end
    end
    fclose(fid);

    fprintf('Saved electrode analysis outputs to %s\n', save_dir);
end
