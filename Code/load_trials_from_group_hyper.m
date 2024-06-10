%% LOAD_TRIALS_FROM_GROUP_HYPER
% This function loads standard and deviant trials from the GroupHyper
% structures, created by Matt's preprocessing of Drosophila LFP data

function [dvt,std] = load_trials_from_group_hyper(overVar_file,...
    deviant_group_numbers, standard_group_numbers, corrected,resample_rate)

% Get the deviant and standard groups
groupHyper = load(overVar_file, 'groupHyper');
for num_dev = 1:length(deviant_group_numbers) % Get all deviant trials
    dvt_group = groupHyper.groupHyper(deviant_group_numbers(num_dev));
    if corrected == 1 %Select filtered data or not
        d_with_nan = dvt_group.DatasCorr;
    elseif corrected == 0
        d_with_nan = dvt_group.Datas;
    end
    if num_dev == 1
        all_d_with_nan = d_with_nan;
    else
        all_d_with_nan = cat(3,all_d_with_nan,d_with_nan);
    end
end
for num_std = 1:length(standard_group_numbers) %Get all standard trials
    std_group = groupHyper.groupHyper(standard_group_numbers(num_std));
    if corrected == 1 %Select filtered data or not
        s_with_nan = std_group.DatasCorr;
    elseif corrected == 0
        s_with_nan = std_group.Datas;
    end
    if num_std == 1
        all_s_with_nan = s_with_nan;
    else
        all_s_with_nan = cat(3,all_s_with_nan,s_with_nan);
    end
end
    



    
% Check if a channel has only NaN values before adding data
dvt.data = all_d_with_nan(all(all(~isnan(all_d_with_nan),3),2),:,:); 
std.data = all_s_with_nan(all(all(~isnan(all_s_with_nan),3),2),:,:);
    

    
dvt.nbchan = size(dvt.data,1);
std.nbchan = size(std.data,1);
dvt.trials = size(dvt.data,3);
std.trials = size(std.data,3);
dvt.times = (1:size(dvt_group.Datas,2))*(1000/resample_rate);
std.times = (1:size(dvt_group.Datas,2))*(1000/resample_rate);
std.chanlocs = [];
dvt.chanlocs = [];

