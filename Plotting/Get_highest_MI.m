%% GET GETS THE ELECTRODE WITH THE HIGHEST MI PEAK FOR TEMPORAL AND FRONTAL AND RETURNS THE ELECTRODE'S DATA IN A STRUCT
function [Highest] = Get_highest_MI(basefold,datatype,condition,partname, cutoff, times)

%Get data in correct format
[tempFFiall, ~, temp_all_mask, ~, ~, ~,tempFFmi1, ~,~,~, ~, ~, ~, ~, ~, ...
          frontFFiall, ~, front_all_mask, ~, ~, ~, frontFFmi1, ~, ~, ~,...
          ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = get_plotting_CoI(basefold,datatype, partname, condition, cutoff, times);

%% FIND ELECTRODE WITH HIGHEST MI FOR TEMPORAL AND FRONTAL
temporal_high_MI_row = max(tempFFmi1);
frontal_high_MI_row = max(frontFFmi1);
temporal_highest = max(temporal_high_MI_row);
frontal_highest = max(frontal_high_MI_row);
temporal_high_MI_row = temporal_high_MI_row == temporal_highest;
frontal_high_MI_row = frontal_high_MI_row == frontal_highest;

%Get all relevant data for the highest electrode
%Temporal
Highest.temp.data = tempFFiall(:,:,temporal_high_MI_row');
Highest.temp.MI = tempFFmi1(:,temporal_high_MI_row);
Highest.temp.sigMask = temp_all_mask(:,:,temporal_high_MI_row');

%Frontal
Highest.front.data = frontFFiall(:,:,frontal_high_MI_row');
Highest.front.MI = frontFFmi1(:,frontal_high_MI_row);
Highest.front.sigMask = front_all_mask(:,:,frontal_high_MI_row');

end






