function CMI_stat = main_CMI_jitterERP(jitterERPFile, outputDir, flyToAnalyseArg, nPermArg)
%% CONDITIONAL MUTUAL INFORMATION FOR JITTER ERP DATA
%
% This is the jitter-ERP equivalent of main_MI_ERP, but instead of asking
% whether the ERP contains information about stimulus class, it asks:
%
%   I(response ; surprise | number_of_preceding_standards)
%
% In words:
%   Does the single-trial neural response contain information about
%   surprise after controlling for the immediate repetition history?
%
% The immediate repetition history, number_of_preceding_standards, is the
% control for stimulus-specific adaptation / run-length effects.
%
% The default surprise variable is recentStimWindowSurpriseBits, which was
% calculated from the previous N valid stimuli. 

%% ---------------- USER SETTINGS ----------------

% Path to Robin Ince's GCMI functions. The script chooses the cluster path
% if it exists, otherwise it uses the local Windows path.
localGcmiPath = 'C:\Users\maxim\OneDrive\Desktop\PhD\Juho Code\Co-I Pipeline\GCMI_master\matlab';
clusterGcmiPath = '/data/SBBS-PIDProject/Maxime/GCMI_master/matlab';
gcmiPath = localGcmiPath;
if exist(clusterGcmiPath,'dir') == 7
    gcmiPath = clusterGcmiPath;
end

% Which trialwise surprise matrix to use from jitterERP.
% Good first choice:
%   'recentStimWindowSurpriseBits'
% Other possible choices:
%   'localWindowSurpriseBits'
%   'localExpSurpriseBits'
surpriseField = 'recentStimWindowSurpriseBits';

% Which column of the surprise matrix to use.
% For recentStimWindowSurpriseBits this refers to recentStimWindowTrials.
% For localWindowSurpriseBits this refers to localWindowTrials.
% For localExpSurpriseBits this refers to localTauTrials.
windowValue = 10;

% Use 'all', 'deviant', or 'standard'.
% 'all' asks whether response tracks surprise across all jittering trials.
% 'deviant' asks the cleaner prediction-error question for deviant events.
trialSelection = 'all';

% If true, the conditioning variable becomes:
%   number_of_preceding_standards + physical_stimulus
% This is stricter because it controls for colour/stimulus identity too.
conditionOnStimulus = false;

% Optional per-trial baseline correction.
% Leave empty to use the ERP snippets as saved.
% Example if samples 1:100 are pre-stimulus:
%   baselineSamples = 1:100;
baselineSamples = 1:5;

% Conditions with fewer than this many trials are discarded before CMI.
minTrialsPerCondition = 5;

% Number of condition-preserving permutations.
% For each permutation, surprise labels are shuffled only within each
% conditioning level, e.g. within each number_of_preceding_standards bin.
% Set to 0 for a quick no-permutation run.
nPerm = 10;
permAlpha = 0.05;

% Analyse one jitterERP entry, or all entries.
% Use [] for all flies/entries.
% Example:
%   flyToAnalyse = 3;
flyToAnalyse = 1;

% Fixed-window surprise has repeated values because it is based on counts
% from a short history window. Treating it as discrete avoids repeated-value
% warnings from gccmi_ccd and better matches the variable we constructed.
%
% Use:
%   'discreteSurprise'   for fixed-window recent-stimulus surprise
%   'continuousSurprise' for gccmi_ccd
estimatorType = 'discreteSurprise';

%% ---------------- DEFAULT PATHS ----------------

if nargin < 1 || isempty(jitterERPFile)
    localJitterERPFile = fullfile(fileparts(mfilename('fullpath')), ...
        '..','Data','Drosophila_LFP', ...
        'Rdaytime_jitteringExperiment_jitteringExperiment_jitterERP.mat');
    clusterJitterERPFile = '/data/SBBS-PIDProject/Maxime/Prediction_Errors/Data/Drosophila_LFP/Rdaytime_jitteringExperiment_jitteringExperiment_jitterERP.mat';
    jitterERPFile = localJitterERPFile;
    if exist(clusterJitterERPFile,'file') == 2
        jitterERPFile = clusterJitterERPFile;
    end
end

if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(fileparts(jitterERPFile),'CMI_Data');
end

if nargin >= 3 && ~isempty(flyToAnalyseArg)
    flyToAnalyse = flyToAnalyseArg;
end

if nargin >= 4 && ~isempty(nPermArg)
    nPerm = nPermArg;
end

mkdir(outputDir);
addpath(gcmiPath);

%% ---------------- LOAD DATA ----------------

loaded = load(jitterERPFile,'jitterERP');
jitterERP = loaded.jitterERP;

CMI_stat = struct;
CMI_stat.sourceFile = jitterERPFile;
CMI_stat.surpriseField = surpriseField;
CMI_stat.windowValue = windowValue;
CMI_stat.trialSelection = trialSelection;
CMI_stat.conditionOnStimulus = conditionOnStimulus;
CMI_stat.baselineSamples = baselineSamples;
CMI_stat.minTrialsPerCondition = minTrialsPerCondition;
CMI_stat.nPerm = nPerm;
CMI_stat.permAlpha = permAlpha;
CMI_stat.flyToAnalyse = flyToAnalyse;
CMI_stat.estimatorType = estimatorType;
CMI_stat.measure = 'I(response ; surprise | number_of_preceding_standards)';

if conditionOnStimulus
    CMI_stat.measure = 'I(response ; surprise | number_of_preceding_standards, physical_stimulus)';
end

%% ---------------- MAIN ANALYSIS LOOP ----------------

if isempty(flyToAnalyse)
    flyIndsToAnalyse = 1:numel(jitterERP);
else
    flyIndsToAnalyse = flyToAnalyse;
end

for flyInd = flyIndsToAnalyse
    disp(['Calculating CMI for jitterERP entry ',num2str(flyInd),' / ',num2str(numel(jitterERP))])

    % Datas is channels x timepoints x trials.
    % CMI is calculated independently for every channel and timepoint.
    data = double(jitterERP(flyInd).Datas);

    % Optional baseline correction, performed trial by trial.
    if ~isempty(baselineSamples)
        baseline = squeeze(mean(data(:,baselineSamples,:),2,'omitnan'));
        data = data - reshape(baseline,size(data,1),1,size(data,3));
    end

    % Pick the chosen surprise column. Each surprise field is trials x
    % parameter-values, so we first find which column corresponds to the
    % requested memory window / tau.
    windowValues = getWindowValues(jitterERP(flyInd),surpriseField);
    windowCol = find(windowValues == windowValue,1);
    surpriseMatrix = jitterERP(flyInd).(surpriseField);
    surprise = surpriseMatrix(:,windowCol);

    % These are the trial-level variables aligned to Datas(:,:,trial).
    nPrecedingStandards = jitterERP(flyInd).trialMeta(:,5);
    trialType = jitterERP(flyInd).trialMeta(:,4);

    % Decide which trials enter the analysis.
    % The condition variable is the thing we control for in the CMI.
    validTrials = ~isnan(surprise) & ~isnan(nPrecedingStandards);
    if strcmpi(trialSelection,'deviant')
        validTrials = validTrials & trialType == 1;
    elseif strcmpi(trialSelection,'standard')
        validTrials = validTrials & trialType == 0;
    end

    conditionLabels = makeConditionLabels( ...
        nPrecedingStandards, ...
        jitterERP(flyInd).physicalStimulus, ...
        conditionOnStimulus);

    % The estimator call is deliberately kept in a separate function.
    % Later, permutation testing can be added there without changing the
    % loading and trial-selection logic in this main script.
    if strcmpi(estimatorType,'continuousSurprise')
        [cmi,pooledMI,keptTrials,conditionLabelsUsed] = compute_CMI_timecourse_gccmi_ccd( ...
            data, ...
            surprise, ...
            conditionLabels, ...
            validTrials, ...
            minTrialsPerCondition);
        cmiPerm = [];
        pValues = [];
        sigMask = [];
    else
        [cmi,pooledMI,keptTrials,conditionLabelsUsed,cmiPerm,pValues,sigMask] = compute_CMI_timecourse_gcmi_model_cd_conditioned( ...
            data, ...
            surprise, ...
            conditionLabels, ...
            validTrials, ...
            minTrialsPerCondition, ...
            nPerm, ...
            permAlpha);
    end

    CMI_stat.fly(flyInd).dataset = jitterERP(flyInd).dataset;
    CMI_stat.fly(flyInd).flyNum = jitterERP(flyInd).flyNum;
    CMI_stat.fly(flyInd).blockNum = jitterERP(flyInd).blockNum;
    CMI_stat.fly(flyInd).CMI = cmi;
    CMI_stat.fly(flyInd).pooledMI = pooledMI;
    CMI_stat.fly(flyInd).CMI_perm = cmiPerm;
    CMI_stat.fly(flyInd).pValues = pValues;
    CMI_stat.fly(flyInd).sigMask = sigMask;
    CMI_stat.fly(flyInd).nTrialsUsed = sum(keptTrials);
    CMI_stat.fly(flyInd).conditionLabels = conditionLabelsUsed;
    CMI_stat.fly(flyInd).surprise = surprise(keptTrials);
    CMI_stat.fly(flyInd).nPrecedingStandards = nPrecedingStandards(keptTrials);
    CMI_stat.fly(flyInd).trialType = trialType(keptTrials);
end

%% ---------------- SAVE ----------------

if isempty(flyToAnalyse)
    flyTag = 'allflies';
else
    flyTag = ['fly',num2str(flyToAnalyse)];
end

saveName = fullfile(outputDir, ...
    ['CMI_',surpriseField,'_',num2str(windowValue),'_',trialSelection,'_',flyTag,'.mat']);

save(saveName,'CMI_stat','-v7.3');
disp(['Saved CMI results to ',saveName])
end

function windowValues = getWindowValues(thisJitterERP,surpriseField)
% Return the parameter vector that defines the columns of the chosen
% surprise matrix.
switch surpriseField
    case 'recentStimWindowSurpriseBits'
        windowValues = thisJitterERP.recentStimWindowTrials;
    case 'localWindowSurpriseBits'
        windowValues = thisJitterERP.localWindowTrials;
    case 'localExpSurpriseBits'
        windowValues = thisJitterERP.localTauTrials;
end
end

function conditionLabels = makeConditionLabels(nPrecedingStandards,physicalStimulus,conditionOnStimulus)
% Build the discrete conditioning labels for gccmi_ccd.
%
% If conditionOnStimulus is false:
%   condition = number_of_preceding_standards
%
% If conditionOnStimulus is true:
%   condition = number_of_preceding_standards + physical_stimulus
%
% The second option asks a stricter question because it controls for
% stimulus identity as well as immediate repetition history.

conditionLabels = string(nPrecedingStandards);

if conditionOnStimulus
    conditionLabels = conditionLabels + "_" + string(physicalStimulus);
end
end
