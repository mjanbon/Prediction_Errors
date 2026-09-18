function analyze_synergy_redundancy_by_electrode(fly_list, basefold, datatype, activity_tag, condition)
% ANALYZE_SYNERGY_REDUNDANCY_BY_ELECTRODE
% For each electrode (E1..En) compute the average CoI matrix across all
% pairs that include that electrode and itself, then compute total and peak synergy
% (negative values) and redundancy (positive values). Aggregate across
% flies, plot the metrics across electrode positions (retina -> central: Emax->E1),
% save figures, CSV and fit statistics.

    fprintf('Analyzing synergy and redundancy by electrode across %d flies\n', length(fly_list));

    % First pass: discover maximum electrode number across flies by
    % reading the electrodes-of-interest (EoI) file — consistent with
    % how `max_get_plotting_CoI_rf.m` builds self-pair names like 'E3_E3'.
    max_elec_num = 0;
    electrode_set = [];

    % Load EoI and extract electrode list exactly as max_get_plotting_CoI_rf.m
    eoi_file = fullfile(basefold, 'DataEOI', ['EoI_data_', datatype, '.mat']);
    if ~exist(eoi_file, 'file')
        error('EoI file not found: %s', eoi_file);
    end
    load(eoi_file, 'EoI');

    % Build electrode set from EoI.(participant).(activity_tag).(condition)
    % Keep this simple and follow max_get_plotting_CoI_rf: expect
    % EoI.(participant).(activity_tag).(condition) to be an ordered list of
    % electrode labels (e.g., {'E1','E2',...}). Print one example for
    % inspection, then extract numeric electrode numbers by simple parsing
    % of the label (no regex).
    % Minimal check: no verbose example printing

    for f = 1:length(fly_list)
        participant_name = fly_list{f};
        if ~isfield(EoI, participant_name)
            continue;
        end
        pentry = EoI.(participant_name);
        if isstruct(pentry) && isfield(pentry, activity_tag) && isfield(pentry.(activity_tag), condition)
            all_elecs = pentry.(activity_tag).(condition);
        elseif isstruct(pentry) && isfield(pentry, condition)
            all_elecs = pentry.(condition);
        else
            continue;
        end

        % Expect all_elecs to be a cell array of strings like {'E1','E2',...}
        for ei = 1:numel(all_elecs)
            el = all_elecs(ei);
            ename = char(el);
            % Simple parsing: if label starts with 'E' drop it and parse the rest
            if startsWith(ename, 'E')
                num = str2double(ename(2:end));
            else
                num = str2double(ename);
            end
            if ~isnan(num) && num>0
                electrode_set(end+1) = num; %#ok<AGROW>
                max_elec_num = max(max_elec_num, num);
            end
        end
    end

    electrode_set = unique(electrode_set);
    if isempty(electrode_set) || max_elec_num == 0
        error('No electrodes found in EoI file for the provided flies/conditions.');
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

            % Discover this fly's available self-pair fields from CoI data
            if ~isfield(CoI, participant_name) || ~isfield(CoI.(participant_name), condition) || ...
               ~isfield(CoI.(participant_name).(condition), 'data')
                fprintf('    No CoI data for %s - skipping\n', participant_name);
                continue;
            end

            % Proceed without verbose diagnostics of data field names

            % Use the electrodes-of-interest (EoI) loaded earlier (same as
            % max_get_plotting_CoI_rf.m) and build the electrode name list.
            if isfield(EoI, participant_name) && isfield(EoI.(participant_name), activity_tag) && isfield(EoI.(participant_name).(activity_tag), condition)
                all_elecs = EoI.(participant_name).(activity_tag).(condition);
            else
                fprintf('    No EoI entry for %s - skipping\n', participant_name);
                continue;
            end

            % Follow the exact procedure used in max_get_plotting_CoI_rf:
            % iterate the EoI-provided electrode list and look for a self-pair
            % field named '<E>_<E>' (e.g. 'E3_E3'). Use myIsField to detect the
            % field and map the electrode label to its numeric index by simple
            % parsing of the label.
            for ei = 1:numel(all_elecs)
                el = all_elecs(ei);
                ename = char(el);

                % Map label to numeric electrode index (E# -> #)
                if startsWith(ename, 'E')
                    enum = str2double(ename(2:end));
                else
                    enum = str2double(ename);
                end
                if isnan(enum)
                    continue;
                end

                accum = [];
                cnt = 0;

                elec_field = [ename, '_', ename];
                % Search data fields for an exact self-pair match (e.g. 'E15_E15')
                data_fields_local = fieldnames(CoI.(participant_name).(condition).data);
                for pf = 1:numel(data_fields_local)
                    fname = data_fields_local{pf};
                    if strcmp(fname, elec_field)
                        mat = CoI.(participant_name).(condition).data.(fname);
                        if isempty(accum)
                            accum = zeros(size(mat));
                        end
                        accum = accum + mat;
                        cnt = cnt + 1;
                        break; % found exact self-pair
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
                end
            end
            % Do not print per-fly diagnostic summaries here

        catch
            % report generic error without relying on a named exception variable
            fprintf('    Error while processing %s - skipping\n', participant_name);
        end
    end

    % Aggregate across flies (omit NaNs)
    mean_total_syn = nanmean(total_syn_by_fly, 1);
    % Compute SEM per electrode column robustly (avoid nanstd signature issues)
    sem_total_syn = nan(1, size(total_syn_by_fly,2));
    for c = 1:size(total_syn_by_fly,2)
        vals = total_syn_by_fly(~isnan(total_syn_by_fly(:,c)), c);
        if isempty(vals)
            sem_total_syn(c) = NaN;
        else
            sem_total_syn(c) = std(vals, 0) / sqrt(numel(vals));
        end
    end
    mean_total_red = nanmean(total_red_by_fly, 1);
    sem_total_red = nan(1, size(total_red_by_fly,2));
    for c = 1:size(total_red_by_fly,2)
        vals = total_red_by_fly(~isnan(total_red_by_fly(:,c)), c);
        if isempty(vals)
            sem_total_red(c) = NaN;
        else
            sem_total_red(c) = std(vals, 0) / sqrt(numel(vals));
        end
    end

    mean_peak_syn = nanmean(peak_syn_by_fly, 1);
    sem_peak_syn = nan(1, size(peak_syn_by_fly,2));
    for c = 1:size(peak_syn_by_fly,2)
        vals = peak_syn_by_fly(~isnan(peak_syn_by_fly(:,c)), c);
        if isempty(vals)
            sem_peak_syn(c) = NaN;
        else
            sem_peak_syn(c) = std(vals, 0) / sqrt(numel(vals));
        end
    end
    mean_peak_red = nanmean(peak_red_by_fly, 1);
    sem_peak_red = nan(1, size(peak_red_by_fly,2));
    for c = 1:size(peak_red_by_fly,2)
        vals = peak_red_by_fly(~isnan(peak_red_by_fly(:,c)), c);
        if isempty(vals)
            sem_peak_red(c) = NaN;
        else
            sem_peak_red(c) = std(vals, 0) / sqrt(numel(vals));
        end
    end

    % Pooled SEM across all fly/electrode observations for figure shading
    pooled_total_syn_sem = nan(1, max_elec_num);
    pooled_total_red_sem = nan(1, max_elec_num);
    pooled_peak_syn_sem = nan(1, max_elec_num);
    pooled_peak_red_sem = nan(1, max_elec_num);
    for c = 1:max_elec_num
        vals = total_syn_by_fly(:,c);
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            pooled_total_syn_sem(c) = std(vals, 0) / sqrt(numel(vals));
        end
        vals = total_red_by_fly(:,c);
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            pooled_total_red_sem(c) = std(vals, 0) / sqrt(numel(vals));
        end
        vals = peak_syn_by_fly(:,c);
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            pooled_peak_syn_sem(c) = std(vals, 0) / sqrt(numel(vals));
        end
        vals = peak_red_by_fly(:,c);
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            pooled_peak_red_sem(c) = std(vals, 0) / sqrt(numel(vals));
        end
    end

    % Suppress verbose diagnostic summary now that electrode coverage is complete

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
    % Align SEM values with the plotted y-values (use same ordering)
    y_syn_sem = pooled_total_syn_sem(elec_x_plot);
    y_red_plot = mean_total_red(elec_x_plot);
    y_red_sem = pooled_total_red_sem(elec_x_plot);
    % Replace NaN SEMs with 0 so shading draws only where data exist
    y_syn_sem(isnan(y_syn_sem)) = 0;
    y_red_sem(isnan(y_red_sem)) = 0;
    % Draw shaded error regions first so lines are always visible on top
    shade_error_region(elec_x, y_syn_plot, y_syn_sem, [0 0.4 0.8]);
    shade_error_region(elec_x, y_red_plot, y_red_sem, [0.8 0.2 0.2]);
    h1 = plot(elec_x, y_syn_plot, 'b-', 'LineWidth', 0.3); %#ok<NASGU>
    h2 = plot(elec_x, y_red_plot, 'r-', 'LineWidth', 0.3); %#ok<NASGU>

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
        % Do not plot linear or exponential fit curves (user requested)
        % Fit models are still computed below to populate fit_stats where possible.
        % Leave mdl_lin_syn, mdl_lin_red creation intact and skip plotting.
        % Exponential models are computed only when valid positive data exist.
        % (No plotting commands here.)
        % Note: mdl_lin_syn and mdl_lin_red already exist above.
        if sum(valid_syn) >= 2
            mdl_exp_syn = fitlm(x_fit(valid_syn), log(y_syn_fit(valid_syn)));
        else
            mdl_exp_syn = [];
        end
        if sum(valid_red) >=2
            mdl_exp_red = fitlm(x_fit(valid_red), log(y_red_fit(valid_red)));
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
    % Peaks: use the same electrode ordering as totals for plotting
    y_ps_plot = -mean_peak_syn(elec_x_plot);
    % Align SEM and red peak values with the plotted ordering
    y_ps_sem = pooled_peak_syn_sem(elec_x_plot);
    y_pr_plot = mean_peak_red(elec_x_plot);
    y_pr_sem = pooled_peak_red_sem(elec_x_plot);
    y_ps_sem(isnan(y_ps_sem)) = 0;
    y_pr_sem(isnan(y_pr_sem)) = 0;

    % Draw shaded error regions first so lines are on top
    shade_error_region(elec_x, y_ps_plot, y_ps_sem, [0 0.4 0.8]);
    shade_error_region(elec_x, y_pr_plot, y_pr_sem, [0.8 0.2 0.2]);
    plot(elec_x, y_ps_plot, 'b-', 'LineWidth', 0.3);
    plot(elec_x, y_pr_plot, 'r-', 'LineWidth', 0.3);
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
    out_fig = fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_Combined_%dFlies.fig', activity_tag, condition, n_flies));
    out_png = fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_Combined_%dFlies.png', activity_tag, condition, n_flies));
    out_svg = fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_Combined_%dFlies.svg', activity_tag, condition, n_flies));
    saveas(fig, out_fig, 'fig');
    try
        % Force OpenGL-backed PNG rendering to preserve alpha transparency
        print(fig, out_png, '-opengl', '-dpng', '-r300');
    catch
        % Fallback for environments where print with -opengl isn't supported
        try
            saveas(fig, out_png);
        catch
            warning('analyze_by_elec:saveFailed', 'Could not save PNG using print or saveas.');
        end
    end
    print(fig, out_svg, '-dsvg', '-painters');
    save(fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_fits_%dFlies.mat', activity_tag, condition, n_flies)), 'fit_stats');

    % --- Mixed-effects model: compare electrode slopes for synergy vs redundancy ---
    try
        fprintf('Fitting mixed-effects model (Value ~ Electrode*MetricType + (1|Fly))...\n');

        % Prepare long-form table
        Y_syn = -total_syn_by_fly; % make positive for plotting
        Y_red = total_red_by_fly;
        Fly = {};
        Electrode = [];
        MetricType = {};
        Value = [];
        for f = 1:n_flies
            for e = 1:max_elec_num
                if isfinite(Y_syn(f,e))
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Electrode(end+1,1) = e; %#ok<SAGROW>
                    MetricType{end+1,1} = 'syn'; %#ok<SAGROW>
                    Value(end+1,1) = Y_syn(f,e); %#ok<SAGROW>
                end
                if isfinite(Y_red(f,e))
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Electrode(end+1,1) = e; %#ok<SAGROW>
                    MetricType{end+1,1} = 'red'; %#ok<SAGROW>
                    Value(end+1,1) = Y_red(f,e); %#ok<SAGROW>
                end
            end
        end

        Electrode_rev = max_elec_num + 1 - Electrode(:);
        T = table(categorical(Fly), Electrode_rev, categorical(MetricType), Value(:), 'VariableNames', {'Fly','Electrode','MetricType','Value'});
        % Ensure baseline for MetricType is 'syn'
        T.MetricType = categorical(T.MetricType, {'syn','red'});

        % Fit model (requires Statistics Toolbox)
        lme = fitlme(T, 'Value ~ Electrode*MetricType + (1|Fly)');

        % Print fixed effects
        disp('Mixed-effects model fixed effects:');
        disp(lme.Coefficients);

        % Extract interaction term (Electrode:MetricType_red)
        termName = 'Electrode:MetricType_red';
        coeffs = lme.Coefficients;
        idx = find(strcmp(coeffs.Name, termName));
        if ~isempty(idx)
            est = coeffs.Estimate(idx);
            se = coeffs.SE(idx);
            tstat = coeffs.tStat(idx);
            pval = coeffs.pValue(idx);
            % 95% CI for the coefficient
            try
                ci = coefCI(lme);
                ci_idx = ci(idx,:);
            catch
                ci_idx = [NaN NaN];
            end
            fprintf('Interaction %s: Estimate=%g, SE=%g, t=%g, p=%g\n', termName, est, se, tstat, pval);
            fprintf('95%% CI: [%g, %g]\n', ci_idx(1), ci_idx(2));
        else
            warning('Interaction term %s not found in model coefficients.', termName);
        end

        % Compute fixed-effect predictions and CIs across electrode axis
        beta = lme.Coefficients.Estimate;
        covB = lme.CoefficientCovariance;
        xe = (1:max_elec_num)';
        % build new design for syn and red
        xe_rev = max_elec_num + 1 - xe;
        X_syn = [ones(numel(xe_rev),1), xe_rev, zeros(numel(xe_rev),1), zeros(numel(xe_rev),1)];
        X_red = [ones(numel(xe_rev),1), xe_rev, ones(numel(xe_rev),1), xe_rev];
        Xall = [X_syn; X_red];
        ypred = Xall * beta;
        varPred = sum((Xall * covB) .* Xall, 2);
        sePred = sqrt(varPred);
        ylow = ypred - 1.96 * sePred;
        yup  = ypred + 1.96 * sePred;

        % Save prediction results
        pred.Electrode = repmat(xe_rev,2,1);
        pred.MetricType = [repmat({'syn'},numel(xe),1); repmat({'red'},numel(xe),1)];
        pred.YPred = ypred;
        pred.YLow = ylow;
        pred.YUp  = yup;
        save(fullfile(save_dir, sprintf('%s_%s_MixedModelPredictions_%dFlies.mat', activity_tag, condition, n_flies)), 'lme', 'T', 'pred');

        % Plot predicted fits with CI
        figm = figure('Position',[200 200 700 400]); hold on;
        % syn (first block)
        xs = xe;
        ys = ypred(1:numel(xe)); ys_low = ylow(1:numel(xe)); ys_up = yup(1:numel(xe));
        fill([xs; flipud(xs)], [ys_low; flipud(ys_up)], [0 0.4 0.8], 'EdgeColor','none', 'FaceAlpha',0.25);
        plot(xs, ys, '-','Color',[0 0.2 0.6],'LineWidth',0.3);
        % red (second block)
        xr = xe;
        yr = ypred(numel(xe)+1:end); yr_low = ylow(numel(xe)+1:end); yr_up = yup(numel(xe)+1:end);
        fill([xr; flipud(xr)], [yr_low; flipud(yr_up)], [0.8 0.2 0.2], 'EdgeColor','none', 'FaceAlpha',0.25);
        plot(xr, yr, '-','Color',[0.6 0.1 0.1],'LineWidth',0.3);
        xlabel('Electrode number (E_{retina} -> E_{central})'); ylabel('Value (bits)');
        legend({'syn CI','syn fit','red CI','red fit'}, 'Location','northeast');
        title(sprintf('Mixed-effects model predictions (%s %s)', activity_tag, condition));
        grid off;

        outm = fullfile(save_dir, sprintf('%s_%s_MixedModelPredictions_%dFlies.png', activity_tag, condition, n_flies));
        outm_svg = fullfile(save_dir, sprintf('%s_%s_MixedModelPredictions_%dFlies.svg', activity_tag, condition, n_flies));
        set(findall(figm, 'Type', 'Line'), 'LineWidth', 0.3);
        set(findall(figm, 'Type', 'Axes'), 'XGrid', 'off', 'YGrid', 'off', 'LineWidth', 0.25);
        try
            print(figm, outm, '-opengl', '-dpng', '-r300');
        catch
            saveas(figm, outm);
        end
        print(figm, outm_svg, '-dsvg', '-painters');
        close(figm);

        % Compute and save fixed-effect slopes for synergy and redundancy
        try
            coeffsTbl = lme.Coefficients;
            names = coeffsTbl.Name;
            idxE = find(strcmp(names, 'Electrode'));
            idxIE = find(strcmp(names, 'Electrode:MetricType_red'));
            betaVec = lme.Coefficients.Estimate;
            covB = lme.CoefficientCovariance;
            if ~isempty(idxE)
                slope_syn = betaVec(idxE);
                var_syn = covB(idxE, idxE);
                se_syn = sqrt(var_syn);
                ci_syn = slope_syn + [-1, 1] * 1.96 * se_syn;
            else
                slope_syn = NaN; se_syn = NaN; ci_syn = [NaN NaN];
            end
            if ~isempty(idxIE) && ~isempty(idxE)
                slope_red = slope_syn + betaVec(idxIE);
                var_red = covB(idxE,idxE) + covB(idxIE,idxIE) + 2 * covB(idxE,idxIE);
                se_red = sqrt(var_red);
                ci_red = slope_red + [-1, 1] * 1.96 * se_red;
            else
                slope_red = NaN; se_red = NaN; ci_red = [NaN NaN];
            end
            fprintf('Fixed-effect slopes (mixed model): syn = %.4g (SE=%.4g, 95%% CI [%.4g, %.4g]); red = %.4g (SE=%.4g, 95%% CI [%.4g, %.4g])\n', slope_syn, se_syn, ci_syn(1), ci_syn(2), slope_red, se_red, ci_red(1), ci_red(2));
            % Append slopes to the MixedModelPredictions MAT file
            try
                save(fullfile(save_dir, sprintf('%s_%s_MixedModelPredictions_%dFlies.mat', activity_tag, condition, n_flies)), 'slope_syn', 'se_syn', 'ci_syn', 'slope_red', 'se_red', 'ci_red', '-append');
            catch
                % ignore save failures
            end
        catch ME
            warning('analyze_by_elec:slopeComputeFailed', 'Failed to compute mixed-model slopes: %s', ME.message);
        end

        fprintf('Saved mixed-model outputs to %s\n', save_dir);

    catch ME
        warning('analyze_by_elec:mixedModelFailed', 'Mixed-effects model failed: %s', ME.message);
    end

    % --- Mixed-effects models for peak synergy/red (raw values) ---
    try
        fprintf('Fitting mixed-effects model for peak values (Value ~ Electrode*MetricType + (1|Fly))...\n');

        Y_syn_peak = -peak_syn_by_fly; % make positive for plotting/analysis
        Y_red_peak = peak_red_by_fly;
        Fly = {};
        Electrode = [];
        MetricType = {};
        Value = [];
        for f = 1:n_flies
            for e = 1:max_elec_num
                if isfinite(Y_syn_peak(f,e))
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Electrode(end+1,1) = e; %#ok<SAGROW>
                    MetricType{end+1,1} = 'syn'; %#ok<SAGROW>
                    Value(end+1,1) = Y_syn_peak(f,e); %#ok<SAGROW>
                end
                if isfinite(Y_red_peak(f,e))
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Electrode(end+1,1) = e; %#ok<SAGROW>
                    MetricType{end+1,1} = 'red'; %#ok<SAGROW>
                    Value(end+1,1) = Y_red_peak(f,e); %#ok<SAGROW>
                end
            end
        end

        Electrode_rev = max_elec_num + 1 - Electrode(:);
        Tpeak = table(categorical(Fly), Electrode_rev, categorical(MetricType), Value(:), 'VariableNames', {'Fly','Electrode','MetricType','Value'});
        Tpeak.MetricType = categorical(Tpeak.MetricType, {'syn','red'});

        lme_peak = fitlme(Tpeak, 'Value ~ Electrode*MetricType + (1|Fly)');
        disp('Mixed-effects model (peaks) fixed effects:');
        disp(lme_peak.Coefficients);

        % Compute and save fixed-effect slopes for peak synergy and redundancy
        try
            coeffsTbl = lme_peak.Coefficients;
            names = coeffsTbl.Name;
            idxE = find(strcmp(names, 'Electrode'));
            idxIE = find(strcmp(names, 'Electrode:MetricType_red'));
            betaVec = lme_peak.Coefficients.Estimate;
            covB = lme_peak.CoefficientCovariance;
            if ~isempty(idxE)
                slope_syn_peak = betaVec(idxE);
                var_syn_peak = covB(idxE, idxE);
                se_syn_peak = sqrt(var_syn_peak);
                ci_syn_peak = slope_syn_peak + [-1, 1] * 1.96 * se_syn_peak;
            else
                slope_syn_peak = NaN; se_syn_peak = NaN; ci_syn_peak = [NaN NaN];
            end
            if ~isempty(idxIE) && ~isempty(idxE)
                slope_red_peak = slope_syn_peak + betaVec(idxIE);
                var_red_peak = covB(idxE,idxE) + covB(idxIE,idxIE) + 2 * covB(idxE,idxIE);
                se_red_peak = sqrt(var_red_peak);
                ci_red_peak = slope_red_peak + [-1, 1] * 1.96 * se_red_peak;
            else
                slope_red_peak = NaN; se_red_peak = NaN; ci_red_peak = [NaN NaN];
            end
            fprintf('Fixed-effect slopes (peak mixed model): syn = %.4g (SE=%.4g, 95%% CI [%.4g, %.4g]); red = %.4g (SE=%.4g, 95%% CI [%.4g, %.4g])\n', slope_syn_peak, se_syn_peak, ci_syn_peak(1), ci_syn_peak(2), slope_red_peak, se_red_peak, ci_red_peak(1), ci_red_peak(2));
            try
                save(fullfile(save_dir, sprintf('%s_%s_PeakMixedModel_%dFlies.mat', activity_tag, condition, n_flies)), 'lme_peak', 'Tpeak', 'slope_syn_peak', 'se_syn_peak', 'ci_syn_peak', 'slope_red_peak', 'se_red_peak', 'ci_red_peak');
            catch
                % ignore save failures
            end
        catch ME
            warning('analyze_by_elec:slopeComputeFailedPeak', 'Failed to compute peak mixed-model slopes: %s', ME.message);
        end

        fprintf('Saved peak mixed-model outputs to %s\n', save_dir);

    catch ME
        warning('analyze_by_elec:peakMixedModelFailed', 'Peak mixed-effects model failed: %s', ME.message);
    end

    % --- Ratio/log-ratio analysis (per-fly normalization) ---
    try
        transform_mode = getenv('COI_TRANSFORM_MODE');
        if isempty(transform_mode)
            transform_mode = 'logratio';
        end
        transform_mode = lower(strtrim(transform_mode));
        if ~ismember(transform_mode, {'ratio', 'logratio'})
            transform_mode = 'logratio';
        end
        fprintf('Computing %s relative to each fly''s own retina-most electrode and fitting mixed-effects model...\n', transform_mode);
        Y_syn = -total_syn_by_fly; % positive synergy
        Y_red = total_red_by_fly;

        ratio_syn = nan(size(Y_syn));
        ratio_red = nan(size(Y_red));
        logratio_syn = nan(size(Y_syn));
        logratio_red = nan(size(Y_red));
        ref_syn = nan(n_flies,1); % store each fly's chosen reference electrode
        ref_red = nan(n_flies,1);

        % Per-fly normalization: each fly uses its own retina-most (highest-index) electrode as baseline
        for f = 1:n_flies
            % For syn: find the highest electrode index with finite data
            ridx_syn = find(isfinite(Y_syn(f,:)), 1, 'last');
            if ~isempty(ridx_syn)
                ref_syn(f) = ridx_syn;
                if abs(Y_syn(f,ridx_syn)) > 0
                    ratio_syn(f,:) = Y_syn(f,:) ./ Y_syn(f,ridx_syn);
                    logratio_syn(f,:) = log(ratio_syn(f,:));
                end
            end
            % For red: find the highest electrode index with finite data
            ridx_red = find(isfinite(Y_red(f,:)), 1, 'last');
            if ~isempty(ridx_red)
                ref_red(f) = ridx_red;
                if abs(Y_red(f,ridx_red)) > 0
                    ratio_red(f,:) = Y_red(f,:) ./ Y_red(f,ridx_red);
                    logratio_red(f,:) = log(ratio_red(f,:));
                end
            end
        end

        % Also compute per-fly normalization for peak values (same approach)
        Y_syn_peak = -peak_syn_by_fly;
        Y_red_peak = peak_red_by_fly;
        ratio_peak_syn = nan(size(Y_syn_peak));
        ratio_peak_red = nan(size(Y_red_peak));
        logratio_peak_syn = nan(size(Y_syn_peak));
        logratio_peak_red = nan(size(Y_red_peak));
        ref_peak_syn = nan(n_flies,1);
        ref_peak_red = nan(n_flies,1);
        for f = 1:n_flies
            ridx_syn = find(isfinite(Y_syn_peak(f,:)), 1, 'last');
            if ~isempty(ridx_syn)
                ref_peak_syn(f) = ridx_syn;
                if abs(Y_syn_peak(f,ridx_syn)) > 0
                    ratio_peak_syn(f,:) = Y_syn_peak(f,:) ./ Y_syn_peak(f,ridx_syn);
                    logratio_peak_syn(f,:) = log(ratio_peak_syn(f,:));
                end
            end
            ridx_red = find(isfinite(Y_red_peak(f,:)), 1, 'last');
            if ~isempty(ridx_red)
                ref_peak_red(f) = ridx_red;
                if abs(Y_red_peak(f,ridx_red)) > 0
                    ratio_peak_red(f,:) = Y_red_peak(f,:) ./ Y_red_peak(f,ridx_red);
                    logratio_peak_red(f,:) = log(ratio_peak_red(f,:));
                end
            end
        end

        % Print reference electrode distribution across flies
        fprintf('Per-fly reference electrodes (syn): ');
        for e = 1:max_elec_num
            fprintf('E%d=%d  ', e, sum(ref_syn==e));
        end
        fprintf('\nPer-fly reference electrodes (red): ');
        for e = 1:max_elec_num
            fprintf('E%d=%d  ', e, sum(ref_red==e));
        end
        fprintf('\n');

        % Choose the transform to model and plot
        if strcmpi(transform_mode, 'ratio')
            trans_syn = ratio_syn;
            trans_red = ratio_red;
            y_label_text = 'Ratio to per-fly reference';
            transform_label = 'Ratio';
        else
            trans_syn = logratio_syn;
            trans_red = logratio_red;
            y_label_text = 'Log-ratio to per-fly reference';
            transform_label = 'Log-ratio';
        end

        % Choose transformed peak arrays as well
        if strcmpi(transform_mode, 'ratio')
            trans_peak_syn = ratio_peak_syn;
            trans_peak_red = ratio_peak_red;
        else
            trans_peak_syn = logratio_peak_syn;
            trans_peak_red = logratio_peak_red;
        end

        % Mean and SEM across flies for plotting
        mean_trans_syn = nanmean(trans_syn, 1);
        sem_trans_syn = nanstd(trans_syn, 0, 1) ./ sqrt(sum(~isnan(trans_syn),1));
        mean_trans_red = nanmean(trans_red, 1);
        sem_trans_red = nanstd(trans_red, 0, 1) ./ sqrt(sum(~isnan(trans_red),1));

        % Detailed per-fly diagnostics for troubleshooting percent calculation
        % Print raw Y_syn/Y_red at each fly's own ref and at a couple of example electrodes (E1,E4)
        echeck = [1, 4];
        fprintf('\nPer-fly raw and %s diagnostics (syn) — per-fly ref, check E%d & E%d:\n', transform_label, echeck(1), echeck(2));
        fprintf('Fly\tRefE\tRefVal\tE1_val\tE1_%s\tE4_val\tE4_%s\n', lower(transform_label), lower(transform_label));
        for f = 1:n_flies
            % syn values
            refE = ref_syn(f);
            if isfinite(refE) && refE>=1 && refE<=size(Y_syn,2)
                refval = Y_syn(f,refE);
            else
                refval = NaN;
            end
            v1 = NaN; t1 = NaN; v4 = NaN; t4 = NaN;
            if isfinite(Y_syn(f,echeck(1)))
                v1 = Y_syn(f,echeck(1)); t1 = trans_syn(f,echeck(1));
            end
            if isfinite(Y_syn(f,echeck(2)))
                v4 = Y_syn(f,echeck(2)); t4 = trans_syn(f,echeck(2));
            end
            fprintf('%s\t%g\t%g\t%g\t%.4f\t%g\t%.4f\n', fly_list{f}, refE, refval, v1, t1, v4, t4);
        end
        fprintf('\nPer-fly raw and %s diagnostics (red) — per-fly ref, check E%d & E%d:\n', transform_label, echeck(1), echeck(2));
        fprintf('Fly\tRefE\tRefVal\tE1_val\tE1_%s\tE4_val\tE4_%s\n', lower(transform_label), lower(transform_label));
        for f = 1:n_flies
            refE = ref_red(f);
            if isfinite(refE) && refE>=1 && refE<=size(Y_red,2)
                refval = Y_red(f,refE);
            else
                refval = NaN;
            end
            v1 = NaN; t1 = NaN; v4 = NaN; t4 = NaN;
            if isfinite(Y_red(f,echeck(1)))
                v1 = Y_red(f,echeck(1)); t1 = trans_red(f,echeck(1));
            end
            if isfinite(Y_red(f,echeck(2)))
                v4 = Y_red(f,echeck(2)); t4 = trans_red(f,echeck(2));
            end
            fprintf('%s\t%g\t%g\t%g\t%.4f\t%g\t%.4f\n', fly_list{f}, refE, refval, v1, t1, v4, t4);
        end

        % Print a short aggregated diagnostic
        fprintf('\n%s diagnostic (first/last 5 indices) — per-fly normalization:\n', transform_label);
        nshow = min(5, max_elec_num);
        for ii = 1:nshow
            idx = ii;
            fprintf('E%d: syn=%.4f (SEM=%.4f), red=%.4f (SEM=%.4f)\n', idx, mean_trans_syn(idx), sem_trans_syn(idx), mean_trans_red(idx), sem_trans_red(idx));
        end
        if max_elec_num > nshow
            fprintf('...\n');
            starti = max(1, max_elec_num - nshow + 1);
            for ii = starti:max_elec_num
                idx = ii;
                fprintf('E%d: syn=%.4f (SEM=%.4f), red=%.4f (SEM=%.4f)\n', idx, mean_trans_syn(idx), sem_trans_syn(idx), mean_trans_red(idx), sem_trans_red(idx));
            end
        end

        % Build long table for transformed values
        Fly = {};
        Electrode = [];
        MetricType = {};
        Value = [];
        for f = 1:n_flies
            for e = 1:max_elec_num
                if isfinite(trans_syn(f,e))
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Electrode(end+1,1) = e; %#ok<SAGROW>
                    MetricType{end+1,1} = 'syn'; %#ok<SAGROW>
                    Value(end+1,1) = trans_syn(f,e); %#ok<SAGROW>
                end
                if isfinite(trans_red(f,e))
                    Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                    Electrode(end+1,1) = e; %#ok<SAGROW>
                    MetricType{end+1,1} = 'red'; %#ok<SAGROW>
                    Value(end+1,1) = trans_red(f,e); %#ok<SAGROW>
                end
            end
        end

        Electrode_rev = max_elec_num + 1 - Electrode(:);
        Tpct = table(categorical(Fly), Electrode_rev, categorical(MetricType), Value(:), 'VariableNames', {'Fly','Electrode','MetricType','Value'});
        Tpct.MetricType = categorical(Tpct.MetricType, {'syn','red'});

        % Fit mixed-effects model on transformed values
        lme_pct = fitlme(Tpct, 'Value ~ Electrode*MetricType + (1|Fly)');
        disp(['Mixed-effects model on ', transform_label, ' values — fixed effects:']);
        disp(lme_pct.Coefficients);

        % Interaction term for transformed model
        termName = 'Electrode:MetricType_red';
        coeffs = lme_pct.Coefficients;
        idx = find(strcmp(coeffs.Name, termName));
        if ~isempty(idx)
            fprintf('%s interaction %s: Estimate=%g, SE=%g, t=%g, p=%g\n', transform_label, termName, coeffs.Estimate(idx), coeffs.SE(idx), coeffs.tStat(idx), coeffs.pValue(idx));
            try
                ci_pct = coefCI(lme_pct);
                fprintf('95%% CI: [%g, %g]\n', ci_pct(idx,1), ci_pct(idx,2));
            catch
                % ignore
            end
        end

        % Save transformed model outputs (include per-fly refs)
        save(fullfile(save_dir, sprintf('%s_%s_%sMixedModel_%dFlies.mat', activity_tag, condition, transform_label, n_flies)), 'lme_pct', 'Tpct', 'mean_trans_syn', 'sem_trans_syn', 'mean_trans_red', 'sem_trans_red', 'ref_syn', 'ref_red', 'transform_mode');

        % Plot transformed-model predictions (fixed-effect only)
        beta = lme_pct.Coefficients.Estimate;
        covB = lme_pct.CoefficientCovariance;
        xe = (1:max_elec_num)';
        xe_rev = max_elec_num + 1 - xe;
        X_syn = [ones(numel(xe_rev),1), xe_rev, zeros(numel(xe_rev),1), zeros(numel(xe_rev),1)];
        X_red = [ones(numel(xe_rev),1), xe_rev, ones(numel(xe_rev),1), xe_rev];
        Xall = [X_syn; X_red];
        ypred = Xall * beta;
        varPred = sum((Xall * covB) .* Xall, 2);
        sePred = sqrt(varPred);
        ylow = ypred - 1.96 * sePred;
        yup  = ypred + 1.96 * sePred;

        figp = figure('Position',[200 200 700 400]); hold on;
        xs = xe;
        ys = ypred(1:numel(xe)); ys_low = ylow(1:numel(xe)); ys_up = yup(1:numel(xe));
        fill([xs; flipud(xs)], [ys_low; flipud(ys_up)], [0 0.4 0.8], 'EdgeColor','none', 'FaceAlpha',0.25);
        plot(xs, ys, '-','Color',[0 0.2 0.6],'LineWidth',0.3);
        xr = xe;
        yr = ypred(numel(xe)+1:end); yr_low = ylow(numel(xe)+1:end); yr_up = yup(numel(xe)+1:end);
        fill([xr; flipud(xr)], [yr_low; flipud(yr_up)], [0.8 0.2 0.2], 'EdgeColor','none', 'FaceAlpha',0.25);
        plot(xr, yr, '-','Color',[0.6 0.1 0.1],'LineWidth',0.3);
        xlabel('Electrode number (E_{retina} -> E_{central})'); ylabel(y_label_text);
        legend({'syn CI','syn fit','red CI','red fit'}, 'Location','northeast');
        title(sprintf('%s (per-fly normalization) — Mixed-model predictions (%s %s)', transform_label, activity_tag, condition));
        grid off;
        outp = fullfile(save_dir, sprintf('%s_%s_%sMixedModelPredictions_%dFlies.png', activity_tag, condition, transform_label, n_flies));
        outp_svg = fullfile(save_dir, sprintf('%s_%s_%sMixedModelPredictions_%dFlies.svg', activity_tag, condition, transform_label, n_flies));
        set(findall(figp, 'Type', 'Line'), 'LineWidth', 0.3);
        set(findall(figp, 'Type', 'Axes'), 'XGrid', 'off', 'YGrid', 'off', 'LineWidth', 0.25);
        try
            print(figp, outp, '-opengl', '-dpng', '-r300');
        catch
            saveas(figp, outp);
        end
        print(figp, outp_svg, '-dsvg', '-painters');
        close(figp);

        % --- Mixed-effects model on transformed peak values ---
        try
            Fly = {};
            Electrode = [];
            MetricType = {};
            Value = [];
            for f = 1:n_flies
                for e = 1:max_elec_num
                    if isfinite(trans_peak_syn(f,e))
                        Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                        Electrode(end+1,1) = e; %#ok<SAGROW>
                        MetricType{end+1,1} = 'syn'; %#ok<SAGROW>
                        Value(end+1,1) = trans_peak_syn(f,e); %#ok<SAGROW>
                    end
                    if isfinite(trans_peak_red(f,e))
                        Fly{end+1,1} = fly_list{f}; %#ok<SAGROW>
                        Electrode(end+1,1) = e; %#ok<SAGROW>
                        MetricType{end+1,1} = 'red'; %#ok<SAGROW>
                        Value(end+1,1) = trans_peak_red(f,e); %#ok<SAGROW>
                    end
                end
            end

            Electrode_rev = max_elec_num + 1 - Electrode(:);
            Tpct_peak = table(categorical(Fly), Electrode_rev, categorical(MetricType), Value(:), 'VariableNames', {'Fly','Electrode','MetricType','Value'});
            Tpct_peak.MetricType = categorical(Tpct_peak.MetricType, {'syn','red'});

            lme_pct_peak = fitlme(Tpct_peak, 'Value ~ Electrode*MetricType + (1|Fly)');
            disp(['Mixed-effects model on ', transform_label, ' peak values — fixed effects:']);
            disp(lme_pct_peak.Coefficients);

            % Compute and save fixed-effect slopes for transformed peak model
            try
                coeffsTbl = lme_pct_peak.Coefficients;
                names = coeffsTbl.Name;
                idxE = find(strcmp(names, 'Electrode'));
                idxIE = find(strcmp(names, 'Electrode:MetricType_red'));
                betaVec = lme_pct_peak.Coefficients.Estimate;
                covB = lme_pct_peak.CoefficientCovariance;
                if ~isempty(idxE)
                    slope_syn_pct_peak = betaVec(idxE);
                    var_syn_pct_peak = covB(idxE, idxE);
                    se_syn_pct_peak = sqrt(var_syn_pct_peak);
                    ci_syn_pct_peak = slope_syn_pct_peak + [-1, 1] * 1.96 * se_syn_pct_peak;
                else
                    slope_syn_pct_peak = NaN; se_syn_pct_peak = NaN; ci_syn_pct_peak = [NaN NaN];
                end
                if ~isempty(idxIE) && ~isempty(idxE)
                    slope_red_pct_peak = slope_syn_pct_peak + betaVec(idxIE);
                    var_red_pct_peak = covB(idxE,idxE) + covB(idxIE,idxIE) + 2 * covB(idxE,idxIE);
                    se_red_pct_peak = sqrt(var_red_pct_peak);
                    ci_red_pct_peak = slope_red_pct_peak + [-1, 1] * 1.96 * se_red_pct_peak;
                else
                    slope_red_pct_peak = NaN; se_red_pct_peak = NaN; ci_red_pct_peak = [NaN NaN];
                end
                fprintf('Fixed-effect slopes (%s peak mixed model): syn = %.4g (SE=%.4g, 95%% CI [%.4g, %.4g]); red = %.4g (SE=%.4g, 95%% CI [%.4g, %.4g])\n', transform_label, slope_syn_pct_peak, se_syn_pct_peak, ci_syn_pct_peak(1), ci_syn_pct_peak(2), slope_red_pct_peak, se_red_pct_peak, ci_red_pct_peak(1), ci_red_pct_peak(2));
                try
                    save(fullfile(save_dir, sprintf('%s_%s_%sMixedModel_%dFlies.mat', activity_tag, condition, transform_label, n_flies)), 'lme_pct_peak', 'Tpct_peak', 'slope_syn_pct_peak', 'se_syn_pct_peak', 'ci_syn_pct_peak', 'slope_red_pct_peak', 'se_red_pct_peak', 'ci_red_pct_peak', '-append');
                catch
                    % ignore save failures
                end
            catch ME
                warning('analyze_by_elec:slopeComputeFailedPctPeak', 'Failed to compute transformed peak mixed-model slopes: %s', ME.message);
            end

        catch ME
            warning('analyze_by_elec:percentPeakModelFailed', 'Transformed peak mixed-effects model failed: %s', ME.message);
        end

        % Per-fly transformed value from (per-fly ref) to electrode 1 (paired tests)
        delta_syn = trans_syn(:,1);
        delta_red = trans_red(:,1);
        ok = isfinite(delta_syn) & isfinite(delta_red);
        n_ok = sum(ok);
        if n_ok > 1
            % Paired t-test on differences (for comparison)
            d_all = delta_syn(ok) - delta_red(ok);
            [~, p_pair] = ttest(d_all);
            % Wilcoxon signed-rank (nonparametric paired test)
            try
                [p_signrank, ~, stats_signrank] = signrank(d_all);
            catch
                p_signrank = NaN; stats_signrank = struct();
            end

            % Hodges-Lehmann estimator (median) and bootstrap CI for median difference
            medDiff = median(d_all);
            try
                nboot = 5000;
                bootmed = bootstrp(nboot, @(x) median(x), d_all);
                ci_med = prctile(bootmed, [2.5 97.5]);
            catch
                ci_med = [NaN NaN];
            end

            fprintf('\nPer-fly %s (per-fly ref)->E1 (n=%d): mean syn=%.4f, mean red=%.4f\n', transform_label, n_ok, nanmean(delta_syn(ok)), nanmean(delta_red(ok)));
            fprintf(' Paired t-test on differences: p=%.4g\n', p_pair);
            if ~isnan(p_signrank)
                fprintf(' Wilcoxon signed-rank: p=%.4g, V=%g\n', p_signrank, stats_signrank.signedrank);
            else
                fprintf(' Wilcoxon signed-rank: failed to compute\n');
            end
            fprintf(' Median difference (syn-red)=%.2f%%, 95%% CI [%.2f, %.2f]\n', medDiff, ci_med(1), ci_med(2));

            % Save paired test outputs appended to the percent-model MAT (if exists)
            pct_matfile = fullfile(save_dir, sprintf('%s_%s_%sMixedModel_%dFlies.mat', activity_tag, condition, transform_label, n_flies));
            try
                if exist(pct_matfile, 'file')
                    save(pct_matfile, 'delta_syn', 'delta_red', 'd_all', 'p_pair', 'p_signrank', 'medDiff', 'ci_med', '-append');
                else
                    save(pct_matfile, 'delta_syn', 'delta_red', 'd_all', 'p_pair', 'p_signrank', 'medDiff', 'ci_med');
                end
            catch
                warning('analyze_by_elec:SavePairedTestsFailed', 'Failed to save paired test results to %s', pct_matfile);
            end
        else
            fprintf('Not enough flies with valid per-fly ref->E1 data for paired percent tests.\n');
        end

        fprintf('Saved %s-model outputs to %s\n', transform_label, save_dir);

        % Compute and save fixed-effect slopes for transformed mixed model (lme_pct)
        try
            coeffsTbl = lme_pct.Coefficients;
            names = coeffsTbl.Name;
            idxE = find(strcmp(names, 'Electrode'));
            idxIE = find(strcmp(names, 'Electrode:MetricType_red'));
            betaVec = lme_pct.Coefficients.Estimate;
            covB = lme_pct.CoefficientCovariance;
            if ~isempty(idxE)
                slope_syn_pct = betaVec(idxE);
                var_syn_pct = covB(idxE, idxE);
                se_syn_pct = sqrt(var_syn_pct);
                ci_syn_pct = slope_syn_pct + [-1, 1] * 1.96 * se_syn_pct;
            else
                slope_syn_pct = NaN; se_syn_pct = NaN; ci_syn_pct = [NaN NaN];
            end
            if ~isempty(idxIE) && ~isempty(idxE)
                slope_red_pct = slope_syn_pct + betaVec(idxIE);
                var_red_pct = covB(idxE,idxE) + covB(idxIE,idxIE) + 2 * covB(idxE,idxIE);
                se_red_pct = sqrt(var_red_pct);
                ci_red_pct = slope_red_pct + [-1, 1] * 1.96 * se_red_pct;
            else
                slope_red_pct = NaN; se_red_pct = NaN; ci_red_pct = [NaN NaN];
            end
            fprintf('Fixed-effect slopes (%s mixed model): syn = %.4g (SE=%.4g, 95%% CI [%.4g, %.4g]); red = %.4g (SE=%.4g, 95%% CI [%.4g, %.4g])\n', transform_label, slope_syn_pct, se_syn_pct, ci_syn_pct(1), ci_syn_pct(2), slope_red_pct, se_red_pct, ci_red_pct(1), ci_red_pct(2));
            try
                save(fullfile(save_dir, sprintf('%s_%s_%sMixedModel_%dFlies.mat', activity_tag, condition, transform_label, n_flies)), 'slope_syn_pct', 'se_syn_pct', 'ci_syn_pct', 'slope_red_pct', 'se_red_pct', 'ci_red_pct', '-append');
            catch
                % ignore save failures
            end
        catch ME
            warning('analyze_by_elec:slopeComputeFailedPct', 'Failed to compute transformed mixed-model slopes: %s', ME.message);
        end

    catch ME
        warning('analyze_by_elec:percentModelFailed', 'Percent-based analysis failed: %s', ME.message);
    end

    % Build CSV rows
    header = {'Fly','Electrode','ElectrodeNumber','TotalSynergy','TotalRedundancy','PeakSynergy','PeakRedundancy','PairCount'};
    % Initialize csv_rows as a proper 2D cell array; first row is the header
    csv_rows = cell(1, numel(header));
    csv_rows(1, :) = header;
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
        row = csv_rows(r, :);
        for c = 1:length(row)
            val = row{c};
            if isempty(val)
                fprintf(fid, ',');
            elseif isnumeric(val)
                if isnan(val)
                    fprintf(fid, ',');
                else
                    fprintf(fid, '%g,', val);
                end
            else
                fprintf(fid, '%s,', val);
            end
        end
        % replace last comma with newline
        fseek(fid, -1, 'cof');
        fprintf(fid, '\n');
    end
    fclose(fid);

    fprintf('Saved electrode analysis outputs to %s\n', save_dir);
end

% Local simple shade helper (matches the working implementation used in
% analyze_synergy_redundancy_across_flies.m). Defining it here ensures this
% file uses the exact same fill behaviour that produces visible shaded bands.
function h = shade_error_region(x, y, err, color)
    % Ensure row vectors and remove NaN positions so fill doesn't get NaNs
    x = x(:)'; y = y(:)'; err = err(:)';
    valid = isfinite(x) & isfinite(y) & isfinite(err);
    if sum(valid) < 2
        h = [];
        return;
    end
    xv = x(valid);
    yv = y(valid);
    ev = err(valid);

    % Create x coordinates for the error region (forward and back)
    x_polygon = [xv, fliplr(xv)];

    % Create y coordinates for the error region (upper and lower bounds)
    y_upper = yv + ev;
    y_lower = yv - ev;
    y_polygon = [y_upper, fliplr(y_lower)];

    % Plot the error region as a filled polygon (attach to current axes)
    h = fill(x_polygon, y_polygon, color, 'Parent', gca);
    set(h, 'EdgeColor', 'none');
    set(h, 'FaceAlpha', 0.45);
    % Prevent the filled polygon from appearing in legends
    try
        set(h, 'HandleVisibility', 'off');
        lh = get(h, 'Annotation');
        if isstruct(lh) || ~isempty(lh)
            set(get(h, 'Annotation').LegendInformation, 'IconDisplayStyle', 'off');
        end
    catch
        % ignore
    end
end
