function [cmi,pooledMI,keptTrials,conditionLabelsUsed,cmiPerm,pValues,sigMask] = compute_CMI_timecourse_gcmi_model_cd_conditioned(data,surprise,conditionLabels,validTrials,minTrialsPerCondition,nPerm,permAlpha)
%% COMPUTE CMI WITH CONTINUOUS RESPONSES AND DISCRETE SURPRISE
%
% This estimates:
%
%   I(response ; surprise | condition)
%
% where:
%   response  = continuous ERP amplitude
%   surprise  = discrete/repeated surprise values
%   condition = discrete control variable, usually number of preceding standards
%
% The observed CMI is:
%
%   sum_z P(z) * I(response ; surprise | condition = z)
%
% Within each condition, I(response ; surprise) is calculated with Robin
% Ince's gcmi_model_cd estimator.
%
% Permutation testing:
%   surprise labels are shuffled only within each condition. This preserves
%   P(surprise | condition), but breaks the trial-by-trial alignment between
%   response and surprise inside each run-length/stimulus bin.

if nargin < 6 || isempty(nPerm)
    nPerm = 0;
end
if nargin < 7 || isempty(permAlpha)
    permAlpha = 0.05;
end

%% ---------------- TRIAL FILTERING ----------------

surprise = surprise(:);
conditionLabels = conditionLabels(:);
validTrials = validTrials(:);

keptTrials = validTrials & ~isnan(surprise) & conditionLabels ~= "";

conditionList = unique(conditionLabels(keptTrials));
for conditionInd = 1:numel(conditionList)
    thisCondition = conditionList(conditionInd);
    thisMask = keptTrials & conditionLabels == thisCondition;
    if sum(thisMask) < minTrialsPerCondition
        keptTrials(thisMask) = false;
    end
end

conditionLabelsUsed = conditionLabels(keptTrials);
surpriseUsed = surprise(keptTrials);
dataUsed = data(:,:,keptTrials);

cmi = nan(size(data,1),size(data,2));
pooledMI = nan(size(data,1),size(data,2));
cmiPerm = [];
pValues = [];
sigMask = [];

if isempty(surpriseUsed)
    return
end

%% ---------------- LABEL CODING ----------------

% Repeated surprise values become class labels 0...nSurpriseLevels-1.
[~,~,surpriseCode] = unique(surpriseUsed);
surpriseCode = surpriseCode - 1;

% gccmi/gcmi discrete variables use compact labels starting at 0.
[~,~,conditionCode] = unique(conditionLabelsUsed);
conditionCode = conditionCode - 1;

conditionLevels = unique(conditionCode);
nTrials = length(surpriseUsed);

%% ---------------- OBSERVED CMI ----------------

[cmi,pooledMI] = calculateDiscreteConditionalMI( ...
    dataUsed, ...
    surpriseCode, ...
    conditionCode, ...
    conditionLevels, ...
    minTrialsPerCondition, ...
    nTrials, ...
    true);

%% ---------------- CONDITION-PRESERVING PERMUTATIONS ----------------

if nPerm == 0
    return
end

cmiPerm = nan(size(cmi,1),size(cmi,2),nPerm,'single');

for permInd = 1:nPerm
    if mod(permInd,10) == 1 || permInd == nPerm
        disp(['  Permutation ',num2str(permInd),' / ',num2str(nPerm)])
    end

    surprisePerm = surpriseCode;

    for conditionInd = 1:numel(conditionLevels)
        thisCondition = conditionLevels(conditionInd);
        conditionTrialInds = find(conditionCode == thisCondition);
        conditionTrialInds = conditionTrialInds(:);

        % Shuffle only among trials with the same conditioning label.
        shuffledInds = conditionTrialInds(randperm(length(conditionTrialInds)));
        surprisePerm(conditionTrialInds) = surprisePerm(shuffledInds);
    end

    cmiPerm(:,:,permInd) = single(calculateDiscreteConditionalMI( ...
        dataUsed, ...
        surprisePerm, ...
        conditionCode, ...
        conditionLevels, ...
        minTrialsPerCondition, ...
        nTrials, ...
        false));
end

% One-sided p-value: how often was shuffled CMI at least as large as observed?
pValues = (1 + sum(double(cmiPerm) >= cmi,3)) ./ (nPerm + 1);
sigMask = pValues < permAlpha;
end

function [cmi,pooledMI] = calculateDiscreteConditionalMI(dataUsed,surpriseCode,conditionCode,conditionLevels,minTrialsPerCondition,nTrials,showProgress)
% Calculate one channels x timepoints CMI matrix for one fixed assignment
% of surprise labels. Used for both observed and shuffled data.

cmi = nan(size(dataUsed,1),size(dataUsed,2));
pooledMI = nan(size(dataUsed,1),size(dataUsed,2));

for ch = 1:size(dataUsed,1)
    if showProgress
        disp(['  Channel ',num2str(ch),' / ',num2str(size(dataUsed,1))])
    end

    for tIdx = 1:size(dataUsed,2)
        response = squeeze(dataUsed(ch,tIdx,:));
        ok = ~isnan(response);

        cmiThisTimepoint = 0;
        pooledResponse = [];
        pooledSurprise = [];

        for conditionInd = 1:numel(conditionLevels)
            thisCondition = conditionLevels(conditionInd);
            thisMask = ok & conditionCode == thisCondition;

            if sum(thisMask) < minTrialsPerCondition
                continue
            end

            [thisResponse,thisSurprise] = prepareDiscreteEstimatorInputs( ...
                response(thisMask), ...
                surpriseCode(thisMask));

            nSurpriseLevels = max(thisSurprise) + 1;
            if nSurpriseLevels < 2
                continue
            end

            conditionWeight = length(thisResponse) / nTrials;
            try
                cmiThisTimepoint = cmiThisTimepoint + conditionWeight * gcmi_model_cd(thisResponse,thisSurprise,nSurpriseLevels);
            catch
                continue
            end

            pooledResponse = [pooledResponse; thisResponse]; %#ok<AGROW>
            pooledSurprise = [pooledSurprise; thisSurprise]; %#ok<AGROW>
        end

        cmi(ch,tIdx) = cmiThisTimepoint;

        % Rough non-conditional comparison, saved for reference only.
        if ~isempty(pooledSurprise) && numel(unique(pooledSurprise)) > 1
            [pooledResponse,pooledSurprise] = prepareDiscreteEstimatorInputs(pooledResponse,pooledSurprise);
            if max(pooledSurprise) >= 1
                try
                    pooledMI(ch,tIdx) = gcmi_model_cd(pooledResponse,pooledSurprise,max(pooledSurprise)+1);
                catch
                    pooledMI(ch,tIdx) = NaN;
                end
            end
        end
    end
end
end

function [responseOut,labelOut] = prepareDiscreteEstimatorInputs(responseIn,labelIn)
% gcmi_model_cd needs at least two observations in every discrete class and
% non-zero response variance. Sparse surprise levels can happen naturally
% after conditioning on number of preceding standards, so remove them here.

minTrialsPerSurpriseLevel = 2;

responseIn = responseIn(:);
labelIn = labelIn(:);
keep = ~isnan(responseIn) & ~isnan(labelIn);

labels = unique(labelIn(keep));
for labelInd = 1:numel(labels)
    labelMask = keep & labelIn == labels(labelInd);
    if sum(labelMask) < minTrialsPerSurpriseLevel
        keep(labelMask) = false;
    end
end

responseOut = responseIn(keep);
labelOut = labelIn(keep);

if numel(responseOut) < 3 || std(responseOut) == 0
    responseOut = [];
    labelOut = 0;
    return
end

[~,~,labelOut] = unique(labelOut);
labelOut = labelOut - 1;
end
