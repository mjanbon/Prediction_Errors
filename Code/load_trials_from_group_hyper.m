%% LOAD_TRIALS_FROM_GROUP_HYPER
% This function loads 

function [dvt,std] = load_trials_from_group_hyper(overVar_file,...
    deviant_group_number, standard_group_number, corrected,resample_rate)

% Get the deviant and standard groups
groupHyper = load(overVar_file, 'groupHyper');
dvt_group = groupHyper.groupHyper(deviant_group_number);
std_group = groupHyper.groupHyper(standard_group_number);

%Select filtered data or not
if corrected == 1
    d_with_nan = dvt_group.DatasCorr;
    s_with_nan = std_group.DatasCorr;
elseif corrected == 0
    d_with_nan = dvt_group.Datas;
    s_with_nan = std_group.Datas;
end
    
% Check if a channel has only NaN values before adding data
dvt.data = d_with_nan(all(all(~isnan(d_with_nan),3),2),:,:); 
std.data = s_with_nan(all(all(~isnan(s_with_nan),3),2),:,:);
    

    
dvt.nbchan = size(dvt.data,1);
std.nbchan = size(std.data,1);
dvt.trials = size(dvt.data,3);
std.trials = size(std.data,3);
dvt.times = (1:size(dvt_group.Datas,2))*(1000/resample_rate);
std.times = (1:size(dvt_group.Datas,2))*(1000/resample_rate);
std.chanlocs = [];
dvt.chanlocs = [];

