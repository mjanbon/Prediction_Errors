function tests = test_sleep_window_MI_scan
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.oldPath = path;
helpers = fileparts(fileparts(mfilename('fullpath')));
toolboxes = fileparts(fileparts(helpers));
addpath(helpers,fullfile(toolboxes,'CNM'),fullfile(toolboxes,'GCMI_master','matlab'));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testMatchedWindowsAndChannelEquivalence(testCase)
rng(9);
data = {randn(15,10,23),randn(15,10,23)};
s = compute_sleep_window_MI_scan(data,1:2,[15 14],10,4,0,0.05);
verifyEqual(testCase,s.trialnum,10);
verifyEqual(testCase,s.sleep_starts,[1 5 9 13 14]);
verifyEqual(testCase,s.window_position([1 end]),[0 1]);
verifyEqual(testCase,size(s.sleep_MI),[2 10 5]);
verifyFalse(testCase,isfield(s,'wake_MI'));
verifyFalse(testCase,isfield(s,'MI_difference'));
one = compute_sleep_window_MI_scan(data,1:2,14,10,4,0,0.05);
verifyEqual(testCase,one.sleep_MI,s.sleep_MI(2,:,:),'AbsTol',1e-12);
end

function testSingleWindow(testCase)
rng(3);
data = {randn(2,10,10),randn(2,10,10)};
s = compute_sleep_window_MI_scan(data,1:2,1:2,[],[],0,0.05);
verifyEqual(testCase,s.sleep_starts,1);
verifyEqual(testCase,s.window_position,0);
verifyEqual(testCase,s.sleep_MI,s.full_sleep_MI,'AbsTol',1e-12);
end

function testGroupUsesFliesNotWindows(testCase)
folder = tempname;
mkdir(folder);
summary = struct('subject','','channel_list',[15 14],'condition','BSLEEP', ...
    'activity_tag','sleep','baseline',1:5, ...
    'group_numbers',{{[3 4],[9 10]}},'corrected',0, ...
    'time_ms',-25:74,'options',struct('windowTrials',[],'windowStep',[]));
for i = 1:3
    summary.subject = sprintf('Synthetic%d',i);
    if i == 1, x = [0 1]; value = 0.1;
    elseif i == 2, x = linspace(0,1,11); value = 0.3;
    else, x = 0; value = 100; end
    summary.window_position = x;
    summary.sleep_mean = repmat(value,2,numel(x));
    save(fullfile(folder,[summary.subject,'_sleep_window_MI_summary.mat']),'summary');
end
g = plot_sleep_window_MI_scan(folder);
verifyEqual(testCase,g.subjects,{'Synthetic1','Synthetic2'});
verifyEqual(testCase,g.mean_sleep,0.2*ones(2,101),'AbsTol',1e-12);
verifyEqual(testCase,g.sem_sleep,0.1*ones(2,101),'AbsTol',1e-12);
verifyEqual(testCase,g.nFlies,2*ones(2,101));
verifyTrue(testCase,isfile(fullfile(folder,'All_flies_BSLEEP_all_electrodes.svg')));
end
