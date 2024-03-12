%% GET GETS THE ELECTRODE WITH THE HIGHEST MI PEAK FOR TEMPORAL AND FRONTAL AND RETURNS THE ELECTRODE'S DATA IN A STRUCT
function [Highest, MI] = Get_highest_MI(basefold,datatype,condition, participantnum,partname, cut, times, participants)

conditions = {'attention' 'counting'};
for participanti = 1:length(participants)
    for conditioni = 1:2
        condition = char(conditions(conditioni));
        participant = char(participants(participanti));

        if strcmp(participant,'EOW') == 1 && strcmp(condition,'counting') == 1
            continue
        end

        
        %Get data in correct format
        [tempFFiall, ~, temp_all_mask, ~, ~, ~,tempFFmi1, ~,~,~, ~, ~, ~, ~, ~, ...
            frontFFiall, ~, front_all_mask, ~, ~, ~, frontFFmi1, ~, ~, ~,...
            ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, temp_elecs, front_elecs] = get_plotting_CoI(basefold,datatype, participanti, participant, condition, cut, times);

        %% FIND ELECTRODE WITH HIGHEST MI FOR TEMPORAL AND FRONTAL
%         temporal_high_MI_row = max(tempFFmi1);
%         frontal_high_MI_row = max(frontFFmi1);
%         temporal_highest = max(temporal_high_MI_row);
%         frontal_highest = max(frontal_high_MI_row);
%         temporal_high_MI_row = temporal_high_MI_row == temporal_highest;
%         frontal_high_MI_row = frontal_high_MI_row == frontal_highest;
%         temporal_highest_elec = temp_elecs(temporal_high_MI_row);
%         frontal_highest_elec = front_elecs(frontal_high_MI_row);
     
        if length(temporal_highest_elec) > 1
            temporal_highest_elec = temporal_highest_elec(1);
            temporal_high_MI_row(10) = 0;
        end

        if length(frontal_highest_elec) > 1
            frontal_highest_elec = frontal_highest_elec(1);
        end


        %Get all relevant data for the highest electrode
        %Temporal
        Highest.(participant).(condition).temp.data = tempFFiall(:,:,temporal_high_MI_row');
        Highest.(participant).(condition).temp.MI = tempFFmi1(:,temporal_high_MI_row);
        Highest.(participant).(condition).temp.sigMask = temp_all_mask(:,:,temporal_high_MI_row');
        Highest.(participant).(condition).temp.label = char(temporal_highest_elec);
       

        %Frontal
        Highest.(participant).(condition).front.data = frontFFiall(:,:,frontal_high_MI_row');
        Highest.(participant).(condition).front.MI = frontFFmi1(:,frontal_high_MI_row);
        Highest.(participant).(condition).front.sigMask = front_all_mask(:,:,frontal_high_MI_row');
        Highest.(participant).(condition).front.label = char(frontal_highest_elec);

        MI.(participant).(condition).temp.MI = tempFFmi1(:,temporal_high_MI_row);
        MI.(participant).(condition).front.MI = frontFFmi1(:,frontal_high_MI_row);

        if isempty(frontFFmi1) == 1
            MI.(participant).(condition).front.MI(1:length(times)) = 0;
        end

        if isempty(tempFFmi1) == 1
            MI.(participant).(condition).temp.MI(1:length(times)) = 0;
        end

    end
end






