function outputDir = run_sleep_window_MI_comparison(USING_HPC, participantIndices, options)
% All-electrode sliding-window MI scan using the groups selected by Max_get_param.
% run_sleep_window_MI_comparison()       : all configured flies, local data.
% run_sleep_window_MI_comparison(0,1)    : first fly locally.
% run_sleep_window_MI_comparison(2,3)    : third fly on Apocrita.
% Options: channels, windowTrials, windowFraction, windowStep, nWindows, nPerm,
% dataRoot, outputDir. windowFraction=0.5 uses half of the selected
% activity's available balanced trials per sliding window.
if nargin < 1 || isempty(USING_HPC), USING_HPC = 0; end
if nargin < 2, participantIndices = []; end
if nargin < 3, options = struct(); end
defaults = struct('channels',15:-1:1,'windowTrials',[],'windowFraction',1,'windowStep',[],'nWindows',[], ...
    'nPerm',0,'dataRoot','','outputDir','');
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(options,names{k}), options.(names{k}) = defaults.(names{k}); end
end
validateattributes(options.windowFraction,{'numeric'},{'scalar','>',0,'<=',1});
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot,'Code'),fullfile(repoRoot,'Helpers'));
if USING_HPC == 2
    toolboxRoot = '/data/SBBS-PIDProject/Maxime';
else
    toolboxRoot = fileparts(repoRoot);
end
addpath(fullfile(toolboxRoot,'CNM'),fullfile(toolboxRoot,'GCMI_master','matlab'));
[basefold,datatype,~,condition,~,participants,~,srate,activityTag, ...
    deviantGroups,standardGroups,corrected,stimOnset,baseline] = Max_get_param(USING_HPC,0);
if ~isempty(options.dataRoot), basefold = options.dataRoot; end
if isempty(participantIndices), participantIndices = 1:numel(participants); end
validateattributes(participantIndices,{'numeric'},{'vector','integer','positive','<=',numel(participants)});
if isempty(options.outputDir)
    outputDir = fullfile(basefold,'MI_Data','sleep_window_scan_all_electrodes', ...
        [activityTag,'_',condition]);
    if ~isempty(options.windowTrials)
        outputDir = fullfile(outputDir,sprintf('trials_%d',options.windowTrials));
    elseif options.windowFraction ~= 1
        outputDir = fullfile(outputDir,['windowFraction_',strrep(num2str(options.windowFraction),'.','p')]);
    end
    if ~isempty(options.windowStep)
        outputDir = fullfile(outputDir,sprintf('step_%d',options.windowStep));
    end
    if ~isempty(options.nWindows)
        outputDir = fullfile(outputDir,sprintf('windows_%d',options.nWindows));
    end
else
    outputDir = options.outputDir;
end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

completed = 0;
for i = participantIndices(:).'
    subject = participants{i};
    source = fullfile(basefold,datatype,[subject,'_',condition,'.mat']);
    if ~exist(source,'file')
        warning('SleepScan:missing','Skipping %s: missing %s',subject,source);
        continue
    end
    fprintf('\n%s: loading %s\n',subject,source);
    % Keep the same balancing and colour-concatenation as
    % load_trials_from_group_hyper. The selected activity comes from
    % Max_get_param, so set activity_tag there before running.
    loaded = load(source,'groupHyper');
    groups = loaded.groupHyper;
    [data{1},data{2}] = balanced_pair(groups,deviantGroups,standardGroups,corrected);
    clear loaded groups
    keep = true(size(data{1},1),1);
    for j = 1:2
        keep = keep & reshape(all(all(isfinite(data{j}),3),2),[],1);
    end
    originalChannels = find(keep).';
    % E1..E15 follow the retained-row convention used in the existing MI
    % figures. Refuse additional channel losses instead of shifting that map.
    if nnz(keep) ~= 15 || any(cellfun(@(x) size(x,2) ~= 100,data))
        warning('SleepScan:invalidData', ...
            'Skipping %s: expected 15 retained electrodes and 100 samples; got %d electrodes and %s samples.', ...
            subject,nnz(keep),mat2str(cellfun(@(x) size(x,2),data)));
        clear data
        continue
    end
    for j = 1:2, data{j} = data{j}(keep,:,:); end
    rng(i-1,'twister');
    windowTrials = options.windowTrials;
    if isempty(windowTrials)
        windowTrials = floor(options.windowFraction * min(cellfun(@(x) size(x,3),data)));
    end
    summary = compute_sleep_window_MI_scan(data,baseline,options.channels, ...
        windowTrials,options.windowStep,options.nPerm,0.05,options.nWindows);
    clear data
    summary.subject = subject;
    summary.condition = condition;
    summary.activity_tag = activityTag;
    summary.source = source;
    summary.originalChannelNumbers = originalChannels(options.channels);
    summary.time_ms = (0:99)*(1000/srate)-stimOnset;
    summary.baseline = baseline;
    summary.corrected = corrected;
    summary.group_numbers = {deviantGroups,standardGroups};
    summary.options = options;
    summary.trial_order = 'colour-concatenated within deviant and standard classes; not recording time';
    save(fullfile(outputDir,[subject,'_',condition,'_sleep_window_MI_summary.mat']),'summary','-v7.3');
    plot_sleep_window_MI_scan(outputDir,subject);
    completed = completed+1;
end
if completed == 0
    error('SleepScan:noData','No requested flies were processed. See missing-file/sample-count messages above.');
end
% Each cluster array task saves only its own figures; aggregate once all
% tasks finish by calling plot_sleep_window_MI_scan on their shared folder.
if numel(participantIndices) > 1
    plot_sleep_window_MI_scan(outputDir);
end
fprintf('\nCompleted %d/%d requested flies. Outputs: %s\n',completed,numel(participantIndices),outputDir);
end

function [dev,std] = balanced_pair(groups,devGroups,stdGroups,corrected)
allGroups = [devGroups,stdGroups];
n = min(arrayfun(@(g) size(groups(g).Datas,3),allGroups));
if corrected, field = 'DatasCorr'; else, field = 'Datas'; end
dev = []; std = [];
for g = devGroups
    dev = cat(3,dev,groups(g).(field)(:,:,1:n));
end
for g = stdGroups
    std = cat(3,std,groups(g).(field)(:,:,1:n));
end
end
