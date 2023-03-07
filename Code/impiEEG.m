
%% IMPORT DATA
% Imports and preprocesses data from standard and deviant trials for the
% specified datatype, condition and participants


function [iEEG1, iEEG2] = impiEEG(participantnum, basefold, datatype, condition, srate, low_cutoff, high_cutoff, filt_order,re_epoch, dev_epochs, std_epochs, epoch_length)


%% Defines the path(s)
addpath('/home/jma201/coi-ieeg')
addpath('/home/jma201/iEEG')
path = basefold; 
datatype = char(datatype);
type = char(condition);

%all the files in specified folder (datatype)
partfilepath = strcat(path,'/',datatype);
filesinfold = dir(partfilepath);                                           

for filei = 1: length(filesinfold)
    filename = cellstr(filesinfold(filei).name);
    allfilenames (1, filei) = filename;                                     
end

allfilenames = allfilenames (3:end);

%Species path + filename to load the data from
STDfilepath = strcat(path,datatype,'/',allfilenames(participantnum),'/',condition,'/',allfilenames(participantnum),'_', type, '_sta.set')
DVTfilepath = strcat(path,datatype,'/',allfilenames(participantnum),'/',condition,'/',allfilenames(participantnum),'_', type, '_dev.set')

%% LOAD DATA
if isfile(STDfilepath) == 1 &&  isfile(DVTfilepath) == 1

    %Load data for deviant trials and re-epoch if re-epoch ==1
    iEEG1 = pop_loadset(DVTfilepath);
    iEEG1 = pop_resample( iEEG1, srate);
    iEEG1 = pop_eegfiltnew(iEEG1, low_cutoff ,high_cutoff,filt_order,0,[],0);
    if re_epoch == 1
        iEEG1 = pop_epoch(iEEG1, dev_epochs, epoch_length, 'newname', 'dvt', 'epochinfo', 'yes');  %LOCAL
    end
    iEEG1 = eeg_checkset( iEEG1 );

    %Load data for standard trials and re-epoch if re-epoch ==1
    iEEG2 = pop_loadset(STDfilepath);
    iEEG2 = pop_resample(iEEG2, srate);
    iEEG2 = pop_eegfiltnew(iEEG2, low_cutoff ,high_cutoff,filt_order,0,[],0);
    if re_epoch == 1
        iEEG2 = pop_epoch(iEEG2, std_epochs, epoch_length, 'newname', 'dvt', 'epochinfo', 'yes');  %LOCAL
    end
    iEEG2 = eeg_checkset(iEEG2);

else
    iEEG1 = [];
    iEEG2 = [];

end
end
