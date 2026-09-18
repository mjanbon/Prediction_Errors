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

    h1 = plot(elec_x, y_syn_plot, 'b-', 'LineWidth', 2); %#ok<NASGU>
    shade_error_region(elec_x, y_syn_plot, y_syn_sem, [0 0.4 0.8]);
    h2 = plot(elec_x, y_red_plot, 'r-', 'LineWidth', 2); %#ok<NASGU>
    shade_error_region(elec_x, y_red_plot, y_red_sem, [0.8 0.2 0.2]);

    xlabel('Electrode (E_{retina} -> E_{central})');
    xticks(elec_x);
    xticklabels(arrayfun(@(n)sprintf('E%d', n), fliplr(elec_x), 'UniformOutput', false));
    xlim([min(elec_x)-0.5, max(elec_x)+0.5]);
    ylabel('Total Information (bits)');
    title(sprintf('%s %s: Total Synergy & Redundancy by Electrode (n=%d flies)', activity_tag, condition, n_flies));
    legend({'Synergy','Redundancy'}, 'Location', 'northeast', 'Box', 'off');

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

    plot(elec_x, y_ps_plot, 'b-', 'LineWidth', 2);
    shade_error_region(elec_x, y_ps_plot, y_ps_sem, [0 0.4 0.8]);
    plot(elec_x, y_pr_plot, 'r-', 'LineWidth', 2);
    shade_error_region(elec_x, y_pr_plot, y_pr_sem, [0.8 0.2 0.2]);
    xlabel('Electrode (E_{retina} -> E_{central})');
    xticks(elec_x);
    xticklabels(arrayfun(@(n)sprintf('E%d', n), fliplr(elec_x), 'UniformOutput', false));
    xlim([min(elec_x)-0.5, max(elec_x)+0.5]);
    ylabel('Peak Value (bits)');
    title('Peak Synergy & Redundancy by Electrode');

    sgtitle(sprintf('%s %s: Synergy/Redundancy by Electrode (n=%d flies)', activity_tag, condition, n_flies));

    % Save outputs: figure, fits and CSV with per-fly per-electrode metrics
    saveas(fig, fullfile(save_dir, sprintf('%s_%s_SynRedByElectrode_Combined_%dFlies.png', activity_tag, condition, n_flies)));
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