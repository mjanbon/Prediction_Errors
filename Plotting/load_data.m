function [elecs, Co_I] = load_data_monkey_JI_XY()

path = strcat('D:\iEEG\Results\Marmo\Ji\XY\results_Ji_XY');
files = dir(strcat('D:\iEEG\Results\Marmo\Ji\XY\results_Ji_XY'));
files = files(3:end);
cd (path)

for filei = 1 : length(files)
    filename = files(filei).name;
    elecs(filei).names = filename;

end

for eleci = 1: length(elecs)
    load(elecs(eleci).names)
    elec_name = elecs(eleci).names;
    elec_name = erase(elec_name, '.mat');
    new_elec_name = erase(elec_name,'_permuted');
    identifier = erase(elec_name,'Ji_');
    identifier = erase(identifier,'permuted_');
    part_name = erase(new_elec_name, identifier);
    part_name = part_name(1:length(part_name)-1);
    new_elec_name = erase(new_elec_name,part_name);
    new_elec_name = new_elec_name(2:length(new_elec_name));

    Co_I.(part_name).XY.data.(new_elec_name) = CoI.(elec_name).CoI.Ji.XY.data.(identifier);
    Co_I.(part_name).XY.sigMask.(new_elec_name) = CoI.(elec_name).CoI.Ji.XY.sigMask.(identifier);
    Co_I.(part_name).XY.MI1.(new_elec_name) = CoI.(elec_name).CoI.Ji.XY.MI1.(identifier);
    Co_I.(part_name).XY.MI2.(new_elec_name) = CoI.(elec_name).CoI.Ji.XY.MI2.(identifier);
    Co_I.(part_name).XY.joint.(new_elec_name) = CoI.(elec_name).CoI.Ji.XY.joint.(identifier);
end
end