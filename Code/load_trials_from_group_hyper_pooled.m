%% LOAD_TRIALS_FROM_GROUP_HYPER_POOLED
% Load standard and deviant trials from groupHyper by pooling matching
% colour/type groups across time windows, then balancing those pooled bins.
%
% This differs from load_trials_from_group_hyper.m, which first clips every
% selected group to the smallest individual group size. Use this loader when
% the requested group numbers are separate time windows that should be
% combined, such as first + last 2 minutes of sleep, while still preserving
% green/blue balance.

function [dvt, std] = load_trials_from_group_hyper_pooled(overVar_file, ...
    deviant_group_numbers, standard_group_numbers, corrected, resample_rate)

groupHyper = load(overVar_file, 'groupHyper');

deviant_bins = make_pooled_balance_bins(deviant_group_numbers);
standard_bins = make_pooled_balance_bins(standard_group_numbers);
all_bins = [deviant_bins, standard_bins];

short_num_trials = inf;
for bin_idx = 1:numel(all_bins)
    pooled_bin = collect_group_trials(groupHyper.groupHyper, all_bins{bin_idx}, corrected);
    short_num_trials = min(short_num_trials, size(pooled_bin, 3));
end

all_d_with_nan = collect_balanced_bins(groupHyper.groupHyper, deviant_bins, corrected, short_num_trials);
all_s_with_nan = collect_balanced_bins(groupHyper.groupHyper, standard_bins, corrected, short_num_trials);

% Keep only channels that are valid in both pooled deviant and pooled standard
% data. This prevents the channel set from differing between classes.
dvt_valid_channels = all(all(~isnan(all_d_with_nan), 3), 2);
std_valid_channels = all(all(~isnan(all_s_with_nan), 3), 2);
valid_channels = dvt_valid_channels & std_valid_channels;

dvt.data = all_d_with_nan(valid_channels, :, :);
std.data = all_s_with_nan(valid_channels, :, :);

dvt.nbchan = size(dvt.data, 1);
std.nbchan = size(std.data, 1);
dvt.trials = size(dvt.data, 3);
std.trials = size(std.data, 3);

first_group = groupHyper.groupHyper(deviant_group_numbers(1));
dvt.times = (1:size(first_group.Datas, 2)) * (1000 / resample_rate);
std.times = dvt.times;
std.chanlocs = [];
dvt.chanlocs = [];

fprintf('Pooled groupHyper trials from %s\n', overVar_file);
fprintf('  Deviant bins %s: %d pooled balanced trials\n', format_bins_for_print(deviant_bins), dvt.trials);
fprintf('  Standard bins %s: %d pooled balanced trials\n', format_bins_for_print(standard_bins), std.trials);
fprintf('  Trials per pooled colour/type bin: %d\n', short_num_trials);
end

function bins = make_pooled_balance_bins(group_numbers)
if mod(numel(group_numbers), 2) ~= 0
    error('Pooled balancing expects an even number of group numbers.');
end

n_types = numel(group_numbers) / 2;
bins = cell(1, n_types);
for type_idx = 1:n_types
    bins{type_idx} = group_numbers(type_idx:n_types:end);
end
end

function all_data = collect_balanced_bins(groupHyper, bins, corrected, n_trials_per_bin)
all_data = [];
for bin_idx = 1:numel(bins)
    pooled_bin = collect_group_trials(groupHyper, bins{bin_idx}, corrected);
    pooled_bin = pooled_bin(:, :, 1:n_trials_per_bin);
    if isempty(all_data)
        all_data = pooled_bin;
    else
        all_data = cat(3, all_data, pooled_bin);
    end
end
end

function pooled_data = collect_group_trials(groupHyper, group_numbers, corrected)
pooled_data = [];
for group_idx = 1:numel(group_numbers)
    this_group = groupHyper(group_numbers(group_idx));
    if corrected == 1
        this_data = this_group.DatasCorr;
    else
        this_data = this_group.Datas;
    end

    if isempty(pooled_data)
        pooled_data = this_data;
    else
        if size(this_data, 1) ~= size(pooled_data, 1) || size(this_data, 2) ~= size(pooled_data, 2)
            error('Group %d has size %dx%d, expected %dx%d.', ...
                group_numbers(group_idx), size(this_data, 1), size(this_data, 2), ...
                size(pooled_data, 1), size(pooled_data, 2));
        end
        pooled_data = cat(3, pooled_data, this_data);
    end
end
end

function text = format_bins_for_print(bins)
parts = cell(1, numel(bins));
for bin_idx = 1:numel(bins)
    parts{bin_idx} = ['[', strtrim(num2str(bins{bin_idx})), ']'];
end
text = strjoin(parts, ' ');
end
