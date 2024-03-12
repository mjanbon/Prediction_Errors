
%% Loads individual CoI-files computed by the HPC to a single structure for convenient handling and plotting
function [elecs, Co_I] = load_data(basefold,datatype,condition,participant_name)

load_path = char(strcat(basefold,'Results\',datatype,'\',participant_name,'\',condition));
files = dir(load_path);
files = files(3:end);
cd (load_path)

for filei = 1 : length(files)
    filename = files(filei).name;
    elecs(filei).names = filename;
end

for eleci = 1: length(elecs)
    load(elecs(eleci).names);
    elec_name = elecs(eleci).names;
    elec_name = erase(elec_name, '.mat');
    new_elec_name = erase(elec_name,'_permuted');
    identifier = erase(elec_name, char(strcat(participant_name,'_')));
    identifier = erase(identifier,'permuted_');
    part_name = erase(new_elec_name, identifier);
    part_name = part_name(1:length(part_name)-1);
    new_elec_name = erase(new_elec_name,part_name);
    new_elec_name = new_elec_name(2:length(new_elec_name));
    eleci
    Co_I.(part_name).(condition).data.(new_elec_name) = CoI.(elec_name).CoI.(participant_name).(condition).data.(identifier)
    Co_I.(part_name).(condition).sigMask.(new_elec_name) = CoI.(elec_name).CoI.(participant_name).(condition).sigMask.(identifier);
    Co_I.(part_name).(condition).MI1.(new_elec_name) = CoI.(elec_name).CoI.(participant_name).(condition).MI1.(identifier);
    Co_I.(part_name).(condition).MI2.(new_elec_name) = CoI.(elec_name).CoI.(participant_name).(condition).MI2.(identifier);
    Co_I.(part_name).(condition).joint.(new_elec_name) = CoI.(elec_name).CoI.(participant_name).(condition).joint.(identifier);
end
end