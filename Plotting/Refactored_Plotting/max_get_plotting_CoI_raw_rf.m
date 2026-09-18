function CoIData = max_get_plotting_CoI_raw_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, times)
% Returns ALL raw CoI, MI, and mask data for every electrode pair (no averaging or summing)

    [~, CoI] = max_load_data(basefold, datatype, condition, activity_tag, participant_name);

    % Load electrodes of interest (EoI)
    eoi_file = fullfile(basefold, 'DataEOI', ['EoI_data_', datatype, '.mat']);
    load(eoi_file, 'EoI');
    all_elecs = EoI.(participant_name).(activity_tag).(condition);

    % Support scalar cutoff (two-region) or two-element vector cutoff
    % (three-region: central, intermediate, peripheral).
    if isscalar(cutoff)
        central_elecs = all_elecs(1:cutoff-1);
        intermediate_elecs = {};
        peripheral_elecs = all_elecs(cutoff:end);
    elseif numel(cutoff) == 2
        c1 = cutoff(1); c2 = cutoff(2);
        if c1 >= c2
            error('cutoff must be increasing: provide [c1 c2] with c1 < c2');
        end
        central_elecs = all_elecs(1:c1-1);
        intermediate_elecs = all_elecs(c1:c2-1);
        peripheral_elecs = all_elecs(c2:end);
    else
        error('cutoff must be a scalar or two-element vector');
    end

    CoIData = struct();

    % Helper for within-region
    % Create within- and cross-region fields dynamically for 2- or 3-region workflows
    regions = {'central','intermediate','peripheral'};
    groups = {central_elecs, intermediate_elecs, peripheral_elecs};

    % within-region
    for r = 1:numel(regions)
        elecs = groups{r};
        if ~isempty(elecs)
            CoIData.(regions{r}) = process_within_region_raw(CoI, participant_name, condition, elecs, times);
        end
    end

    % cross-region and within-group pair fields
    for i = 1:numel(regions)
        for j = i:numel(regions)
            elecsA = groups{i}; elecsB = groups{j};
            if isempty(elecsA) || isempty(elecsB)
                continue;
            end
            fieldname = sprintf('%s_%s', regions{i}, regions{j});
            if i == j
                CoIData.(fieldname) = process_cross_region_raw(CoI, participant_name, condition, elecsA, elecsB, times, true);
            else
                CoIData.(fieldname) = process_cross_region_raw(CoI, participant_name, condition, elecsA, elecsB, times);
            end
        end
    end

    CoIData.central_elecs = central_elecs;
    CoIData.peripheral_elecs = peripheral_elecs;
end

function regionData = process_within_region_raw(CoI, pname, cond, elecs, times)
    n_elecs = numel(elecs);
    nT = length(times);
    regionData.FFi_all = [];
    regionData.mask = [];
    regionData.mask_r = [];
    regionData.mask_s = [];
    regionData.mi1 = [];
    regionData.mi2 = [];
    pair_count = 0;
    for i = 1:n_elecs
        elec_name = [char(elecs(i)), '_', char(elecs(i))];
        if myIsField(CoI.(pname).(cond).data, elec_name)
            pair_count = pair_count + 1;
            regionData.FFi_all(:,:,pair_count) = CoI.(pname).(cond).data.(elec_name);
            regionData.mask(:,:,pair_count) = logical(CoI.(pname).(cond).sigMask.(elec_name));
            regionData.mask_r(:,:,pair_count) = (CoI.(pname).(cond).sigMask.(elec_name) > 0);
            regionData.mask_s(:,:,pair_count) = (CoI.(pname).(cond).sigMask.(elec_name) < 0);
            
            % MI1
            raw_mi1 = CoI.(pname).(cond).MI1.(elec_name);
            if isrow(raw_mi1)
                raw_mi1 = raw_mi1';
            end
            if size(raw_mi1,2) > 1
                warning('MI1 for %s is not a vector, taking only the first column.', elec_name);
                raw_mi1 = raw_mi1(:,1);
            end
            regionData.mi1(:,pair_count) = raw_mi1;

            % MI2
            raw_mi2 = CoI.(pname).(cond).MI2.(elec_name);
            if isrow(raw_mi2)
                raw_mi2 = raw_mi2';
            end
            if size(raw_mi2,2) > 1
                finite_counts = sum(~isnan(raw_mi2), 1);
                [max_count, col_idx] = max(finite_counts);
                if max_count > 0
                    if sum(finite_counts == max_count) > 1
                        warning('MI2 for %s has multiple populated columns; taking the first with the most finite values.', elec_name);
                    end
                    raw_mi2 = raw_mi2(:, col_idx);
                else
                    warning('MI2 for %s is not a vector, taking only the first column.', elec_name);
                    raw_mi2 = raw_mi2(:,1);
                end
            end
            regionData.mi2(:,pair_count) = raw_mi2;
        end
    end
    regionData.n_pairs = pair_count;
end

function regionData = process_cross_region_raw(CoI, pname, cond, elecsA, elecsB, times, skip_duplicates)
    if nargin < 7
        skip_duplicates = false;
    end
    nT = length(times);
    regionData.FFi_all = [];
    regionData.mask = [];
    regionData.mask_r = [];
    regionData.mask_s = [];
    regionData.mi1 = [];
    regionData.mi2 = [];
    elec_pairs = {};
    pair_count = 0;
    for i = 1:numel(elecsA)
        for j = 1:numel(elecsB)
            E1 = char(elecsA(i)); E2 = char(elecsB(j));
            if strcmp(E1, E2)
                continue;
            end
            elec_name = [E1, '_', E2];
            elec_inv = [E2, '_', E1];
            if skip_duplicates && ismember(elec_inv, elec_pairs)
                continue;
            end
            elec_pairs{end+1} = elec_name;
            if myIsField(CoI.(pname).(cond).data, elec_name)
                pair_count = pair_count + 1;
                regionData.FFi_all(:,:,pair_count) = CoI.(pname).(cond).data.(elec_name);
                regionData.mask(:,:,pair_count) = logical(CoI.(pname).(cond).sigMask.(elec_name));
                regionData.mask_r(:,:,pair_count) = (CoI.(pname).(cond).sigMask.(elec_name) > 0);
                regionData.mask_s(:,:,pair_count) = (CoI.(pname).(cond).sigMask.(elec_name) < 0);
                
                % MI1
                raw_mi1 = CoI.(pname).(cond).MI1.(elec_name);
                if isrow(raw_mi1)
                    raw_mi1 = raw_mi1';
                end
                if size(raw_mi1,2) > 1
                    warning('MI1 for %s is not a vector, taking only the first column.', elec_name);
                    raw_mi1 = raw_mi1(:,1);
                end
                regionData.mi1(:,pair_count) = raw_mi1;

                % MI2
                raw_mi2 = CoI.(pname).(cond).MI2.(elec_name);
                if isrow(raw_mi2)
                    raw_mi2 = raw_mi2';
                end
                if size(raw_mi2,2) > 1
                    finite_counts = sum(~isnan(raw_mi2), 1);
                    [max_count, col_idx] = max(finite_counts);
                    if max_count > 0
                        if sum(finite_counts == max_count) > 1
                            warning('MI2 for %s has multiple populated columns; taking the first with the most finite values.', elec_name);
                        end
                        raw_mi2 = raw_mi2(:, col_idx);
                    else
                        warning('MI2 for %s is not a vector, taking only the first column.', elec_name);
                        raw_mi2 = raw_mi2(:,1);
                    end
                end
                regionData.mi2(:,pair_count) = raw_mi2;
            end
        end
    end
    regionData.n_pairs = pair_count;
end