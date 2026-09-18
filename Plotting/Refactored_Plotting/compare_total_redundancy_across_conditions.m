function stats = compare_total_redundancy_across_conditions(participants, basefold, datatype, cutoff, times, activity_tags)
%COMPARE_TOTAL_REDUNDANCY_ACROSS_CONDITIONS Backward-compatible wrapper.
% Use compare_total_coi_across_conditions for the unified implementation.

if nargin < 6
	activity_tags = [];
end

stats_all = compare_total_coi_across_conditions(participants, basefold, datatype, cutoff, times, activity_tags);
stats = stats_all.redundancy;
end