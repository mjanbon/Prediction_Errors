function summary = compute_sleep_window_MI_scan(data, baseline, channels, windowTrials, windowStep, nPerm, alpha)
% DATA contains {deviant, standard}, each channels x samples x trials, with
% the loader's colour-concatenated order. This is a trial-index sensitivity
% scan, not a chronological sleep-bout scan.
if nargin < 6 || isempty(nPerm), nPerm = 0; end
if nargin < 7 || isempty(alpha), alpha = 0.05; end %#ok<NASGU>
if nPerm > 0
    warning('SleepScan:noPermutations', ...
        'nPerm is ignored: compute_sleep_window_MI_scan now calculates sleep/activity MI only, with no wake difference test.');
end
counts = cellfun(@(x) size(x,3),data);
nTime = size(data{1},2);
if isempty(windowTrials), windowTrials = min(counts); end
validateattributes(windowTrials,{'numeric'},{'scalar','integer','>=',2,'<=',min(counts)});
if isempty(windowStep), windowStep = max(1,floor(windowTrials/10)); end
validateattributes(windowStep,{'numeric'},{'scalar','integer','positive'});
validateattributes(nPerm,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(channels,{'numeric'},{'vector','integer','positive','<=',size(data{1},1)});
if any(cellfun(@(x) size(x,2) ~= nTime || size(x,1) ~= size(data{1},1),data))
    error('SleepScan:shape','Deviant and standard inputs must share electrode and time dimensions.');
end
lastStart = min(counts) - windowTrials + 1;
starts = unique([1:windowStep:lastStart,lastStart]);
nCh = numel(channels);
nWindows = numel(starts);
summary.channel_list = channels;
summary.trialnum = windowTrials;
summary.available_trial_counts = counts;
summary.sleep_starts = starts;
summary.window_centers = starts + floor((windowTrials-1)/2);
summary.sleep_window_step = windowStep;
summary.window_position = zeros(size(starts));
if lastStart > 1
    summary.window_position = (starts-1)/(lastStart-1);
end
summary.sleep_MI = nan(nCh,nTime,nWindows);
summary.full_sleep_MI = nan(nCh,nTime);
stim = [zeros(windowTrials,1);ones(windowTrials,1)];

for e = 1:nCh
    ch = channels(e);
    fprintf('  E%d (%d/%d), %d windows\n',ch,e,nCh,nWindows);
    traces = cell(1,2);
    for j = 1:2
        traces{j} = double(reshape(data{j}(ch,:,:),nTime,[]).');
        if ~isempty(baseline)
            traces{j} = traces{j} - mean(traces{j}(:,baseline),2);
        end
    end
    fullSleep = [traces{1};traces{2}];
    fullLabels = [zeros(size(traces{1},1),1);ones(size(traces{2},1),1)];
    summary.full_sleep_MI(e,:) = cnm_MI_stimtime(fullSleep,fullLabels,0).';
    for w = 1:nWindows
        inds = starts(w):(starts(w)+windowTrials-1);
        sleep = [traces{1}(inds,:);traces{2}(inds,:)];
        sleepMI = cnm_MI_stimtime(sleep,stim,0);
        summary.sleep_MI(e,:,w) = sleepMI(:).';
    end
end
summary.sleep_mean = reshape(mean(summary.sleep_MI,2,'omitnan'),nCh,nWindows);
end
