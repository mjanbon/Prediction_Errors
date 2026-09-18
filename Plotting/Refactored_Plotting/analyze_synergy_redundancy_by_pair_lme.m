function res = analyze_synergy_redundancy_by_pair_lme(fly_list, basefold, datatype, activity_tag, condition, outdir)
% ANALYZE_SYNERGY_REDUNDANCY_BY_PAIR_LME
% Collects per-electrode-pair metrics (total synergy, total redundancy,
% peak synergy, peak redundancy) across all flies and fits linear mixed
% effects models using `fitlme` with random intercepts for `Fly`.
%
% Usage:
%   res = analyze_synergy_redundancy_by_pair_lme(fly_list, basefold, datatype, activity_tag, condition)
%   res = analyze_synergy_redundancy_by_pair_lme(..., outdir)
%
% Outputs:
%   - saves a CSV of the per-pair table and a MAT file with `res` containing
%     the table and fitted models for each metric.

if nargin < 6 || isempty(outdir)
    outdir = fullfile(basefold, 'Results', 'CoI_distance', 'By_pair');
end
if ~exist(outdir, 'dir'), mkdir(outdir); end

pairs = {};
flies = {};
E1 = [];
E2 = [];
TotalSyn = [];
TotalRed = [];
PeakSyn = [];
PeakRed = [];

rowi = 0;
for f = 1:numel(fly_list)
    participant = fly_list{f};
    try
        [~, CoI] = max_load_data(basefold, datatype, condition, activity_tag, participant);
    catch
        warning('Could not load CoI for %s - skipping', participant);
        continue;
    end

    if ~isfield(CoI, participant) || ~isfield(CoI.(participant), condition) || ~isfield(CoI.(participant).(condition), 'data')
        warning('No CoI data found for %s - skipping', participant);
        continue;
    end

    data_fields = fieldnames(CoI.(participant).(condition).data);
    for k = 1:numel(data_fields)
        fname = data_fields{k};
        mat = CoI.(participant).(condition).data.(fname);
        if isempty(mat) || ~isnumeric(mat)
            continue;
        end

        % Parse pair label (expect format 'E#_E#')
        parts = strsplit(fname, '_');
        if numel(parts) >= 2
            p1 = parts{1}; p2 = parts{2};
        else
            p1 = fname; p2 = '';
        end
        % extract numeric electrode ids when possible
        e1 = nan; e2 = nan;
        if ~isempty(p1)
            if startsWith(p1, 'E')
                e1 = str2double(p1(2:end));
            else
                e1 = str2double(p1);
            end
        end
        if ~isempty(p2)
            if startsWith(p2, 'E')
                e2 = str2double(p2(2:end));
            else
                e2 = str2double(p2);
            end
        end

        % compute metrics
        total_syn = sum(mat(mat < 0), 'all');
        total_red = sum(mat(mat > 0), 'all');
        peak_syn = min(mat(:));
        peak_red = max(mat(:));

        rowi = rowi + 1;
        pairs{rowi,1} = fname; %#ok<AGROW>
        flies{rowi,1} = participant; %#ok<AGROW>
        E1(rowi,1) = e1; %#ok<AGROW>
        E2(rowi,1) = e2; %#ok<AGROW>
        TotalSyn(rowi,1) = total_syn; %#ok<AGROW>
        TotalRed(rowi,1) = total_red; %#ok<AGROW>
        PeakSyn(rowi,1) = peak_syn; %#ok<AGROW>
        PeakRed(rowi,1) = peak_red; %#ok<AGROW>
    end
end

% Build table with all metrics as columns
T = table(pairs, flies, E1, E2, TotalSyn, TotalRed, PeakSyn, PeakRed, ...
    'VariableNames', {'Pair', 'Fly', 'E1', 'E2', 'TotalSyn', 'TotalRed', 'PeakSyn', 'PeakRed'});

% Convert Fly to categorical for fitlme
T.Fly = categorical(T.Fly);

% Save the table CSV
csvfile = fullfile(outdir, sprintf('%s_%s_per_pair_metrics.csv', activity_tag, condition));
try
    writetable(T, csvfile);
catch
    warning('Could not write CSV to %s', csvfile);
end

% Fit LME models for each metric using Electrode numbers as fixed effects
metrics = {'TotalSyn', 'TotalRed', 'PeakSyn', 'PeakRed'};
res = struct();
res.table = T;
res.fits = struct();

for m = 1:numel(metrics)
    mname = metrics{m};
    try
        formula = sprintf('%s ~ E1 + E2 + (1|Fly)', mname);
        mdl = fitlme(T, formula);
        res.fits.(mname) = mdl;
    catch ME
        warning('fitlme failed for %s: %s', mname, ME.message);
        res.fits.(mname) = [];
    end
end

% Save results
matfile = fullfile(outdir, sprintf('%s_%s_per_pair_metrics_and_fits.mat', activity_tag, condition));
try
    save(matfile, 'res', '-v7.3');
catch
    warning('Could not save MAT to %s', matfile);
end

fprintf('Per-pair table saved to %s and fits saved to %s\n', csvfile, matfile);

end
