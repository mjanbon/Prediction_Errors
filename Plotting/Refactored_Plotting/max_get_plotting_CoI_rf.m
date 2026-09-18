function CoIData = max_get_plotting_CoI_rf(basefold, datatype, participant_name, activity_tag, condition, cutoff, times, comparison_folder)
% MAX_GET_PLOTTING_COI
% Extracts and organizes co-information (CoI) and mutual information (MI) data
% into structured outputs for plotting.
%
% INPUTS:
%   basefold         - Base folder path containing CoI and EoI data
%   datatype         - String identifying the data type (e.g. 'Marmo_EcoG', 'Human_ECoG')
%   participant_name - Name of the participant (e.g. 'S01')
%   activity_tag     - Subcategory label within participant (e.g. 'Resting', 'Task')
%   condition        - Experimental condition label (e.g. 'XY', 'AB')
%   cutoff           - Integer index separating central and peripheral electrodes
%                      (e.g. cutoff = 10 → 1–9 are central, 10–end are peripheral)
%   times            - Vector of timepoints used in the data (e.g. -25:1:74)
%
% OUTPUT:
%   CoIData          - A structured output containing all processed data fields:
%                      .central, .peripheral, .central_peripheral, .central_central,
%                      .peripheral_peripheral
%                      Each subfield includes:
%                          .mean         - average CoI matrix
%                          .FFi_all      - 3D matrix of all electrode-level CoIs
%                          .mask         - combined significance mask
%                          .mask_r/s     - separate redundancy/synergy masks
%                          .mi1, .mi2    - mutual information vectors
%                          .redundant    - only redundant portions of CoI
%                          .synergetic   - only synergetic portions of CoI
%                          .n            - number of electrodes or pairs included
%                      .central_elecs, .peripheral_elecs - lists of electrode labels

    if nargin < 8 || isempty(comparison_folder)
        comparison_folder = 'Comparisons';
    end

    % Load participant-specific CoI struct
    [~, CoI] = max_load_data(basefold, datatype, condition, activity_tag, participant_name, comparison_folder);

    % Load electrodes of interest (EoI)
    eoi_file = fullfile(basefold, 'DataEOI', ['EoI_data_', datatype, '.mat']);
    load(eoi_file, 'EoI');
    all_elecs = EoI.(participant_name).(activity_tag).(condition);

    % Split electrodes into regions based on cutoff.
    % Backwards compatible: if scalar `cutoff` provided, create two regions
    % (central, peripheral) as before. If a two-element vector is provided
    % e.g. [c1 c2], create three regions: central (1:c1-1), intermediate
    % (c1:c2-1) and peripheral (c2:end).
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

    % Initialize structured output
    CoIData = struct();

    % Process within- and cross-regions. Support 2-region and 3-region
    % workflows by building groups dynamically.
    regions = {'central','intermediate','peripheral'};
    groups = {central_elecs, intermediate_elecs, peripheral_elecs};

    % Create within-region fields for non-empty groups
    for r = 1:numel(regions)
        elecs = groups{r};
        if ~isempty(elecs)
            CoIData.(regions{r}) = process_within_region(CoI, participant_name, condition, elecs, times);
        end
    end

    % Create cross-region and within-region-pair fields
    for i = 1:numel(regions)
        for j = i:numel(regions)
            elecsA = groups{i}; elecsB = groups{j};
            if isempty(elecsA) || isempty(elecsB)
                continue;
            end
            fieldname = sprintf('%s_%s', regions{i}, regions{j});
            if i == j
                % within-group (skip duplicate pairs)
                CoIData.(fieldname) = process_cross_region(CoI, participant_name, condition, elecsA, elecsB, times, true);
            else
                CoIData.(fieldname) = process_cross_region(CoI, participant_name, condition, elecsA, elecsB, times);
            end
        end
    end

    % Store electrode lists
    CoIData.central_elecs = central_elecs;
    CoIData.peripheral_elecs = peripheral_elecs;
end


function regionData = process_within_region(CoI, pname, cond, elecs, times)
    regionData = initialize_region_struct(times);
    n_elecs = numel(elecs);
    for i = 1:n_elecs
        elec_name = [char(elecs(i)), '_', char(elecs(i))];
        if myIsField(CoI.(pname).(cond).data, elec_name)
            regionData.FFi_all = cat(3, regionData.FFi_all, CoI.(pname).(cond).data.(elec_name));
            regionData.mask = regionData.mask + logical(CoI.(pname).(cond).sigMask.(elec_name));
            regionData.mask_r = regionData.mask_r + (CoI.(pname).(cond).sigMask.(elec_name) > 0);
            regionData.mask_s = regionData.mask_s + (CoI.(pname).(cond).sigMask.(elec_name) < 0);
            regionData.mi1 = [regionData.mi1, CoI.(pname).(cond).MI1.(elec_name)];
            regionData.mi2 = [regionData.mi2, CoI.(pname).(cond).MI2.(elec_name)];
        end
    end
    regionData.mean = mean(regionData.FFi_all, 3);
    regionData.redundant = regionData.mean .* (regionData.mean > 0);
    regionData.synergetic = regionData.mean .* (regionData.mean < 0);
    regionData.n = n_elecs;
end

function regionData = process_cross_region(CoI, pname, cond, elecsA, elecsB, times, skip_duplicates)
    if nargin < 7
        skip_duplicates = false;
    end
    regionData = initialize_region_struct(times);
    elec_pairs = {};
    nb = 0;
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
                regionData.FFi_all = cat(3, regionData.FFi_all, CoI.(pname).(cond).data.(elec_name));
                regionData.mask = regionData.mask + logical(CoI.(pname).(cond).sigMask.(elec_name));
                regionData.mask_r = regionData.mask_r + (CoI.(pname).(cond).sigMask.(elec_name) > 0);
                regionData.mask_s = regionData.mask_s + (CoI.(pname).(cond).sigMask.(elec_name) < 0);
                regionData.mi1 = [regionData.mi1, CoI.(pname).(cond).MI1.(elec_name)];
                regionData.mi2 = [regionData.mi2, CoI.(pname).(cond).MI2.(elec_name)];
                nb = nb + 1;
            end
        end
    end
    regionData.mean = mean(regionData.FFi_all, 3);
    regionData.redundant = regionData.mean .* (regionData.mean > 0);
    regionData.synergetic = regionData.mean .* (regionData.mean < 0);
    regionData.n = nb;
end

function S = initialize_region_struct(times)
    n = length(times);
    S = struct('FFi_all', [], ...
               'mask', zeros(n), ...
               'mask_r', zeros(n), ...
               'mask_s', zeros(n), ...
               'mi1', [], ...
               'mi2', [], ...
               'mean', [], ...
               'redundant', [], ...
               'synergetic', [], ...
               'n', 0);
end
