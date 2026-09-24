function group = plot_sleep_window_MI_scan(outputDir, subject)
% Plot saved all-electrode scans; no LFP reload or MI recomputation.
% With SUBJECT, plot that fly only. With just OUTPUTDIR, plot an equally
% weighted across-fly mean +/- SEM on a common window-position grid.
if nargin < 2, subject = ''; end
files = dir(fullfile(outputDir,[subject,'*_sleep_window_MI_summary.mat']));
if isempty(files), error('SleepScan:noSummaries','No scan summaries in %s.',outputDir); end
summaries = cell(1,numel(files));
for i = 1:numel(files)
    loaded = load(fullfile(files(i).folder,files(i).name),'summary');
    summaries{i} = loaded.summary;
    if ~isfield(summaries{i}.options,'windowFraction')
        summaries{i}.options.windowFraction = 1;
    end
    if ~isfield(summaries{i}.options,'nWindows')
        summaries{i}.options.nWindows = [];
    end
end
first = summaries{1};
group = [];
if ~isempty(subject)
    for i = 1:numel(summaries)
        s = summaries{i};
        activityTag = get_activity_tag(s);
        draw_panels(s.window_centers,s.sleep_mean,[],get_random_sample_mean(s), ...
            s.channel_list,sprintf('%s: %s, N = %d trials per class', ...
            s.subject,activityTag,s.trialnum), ...
            'Window centre in colour-concatenated selected-activity trials', ...
            fullfile(outputDir,[s.subject,'_',s.condition,'_all_electrodes']));
    end
    return
end
grid = linspace(0,1,101);
channels = first.channel_list;
sleep = nan(numel(channels),numel(grid),numel(files));
randomSample = nan(numel(channels),numel(files));
included = false(1,numel(files));
for i = 1:numel(summaries)
    s = summaries{i};
    if ~isequal(s.channel_list,channels) || ~isequal(s.baseline,first.baseline) || ...
            ~isequal(s.group_numbers,first.group_numbers) || s.corrected ~= first.corrected || ...
            ~isequal(s.time_ms,first.time_ms) || ...
            ~isequal(s.options.windowTrials,first.options.windowTrials) || ...
            s.options.windowFraction ~= first.options.windowFraction || ...
            ~isequal(s.options.windowStep,first.options.windowStep) || ...
            ~isequal(s.options.nWindows,first.options.nWindows)
        error('SleepScan:mixedSettings','Incompatible scan settings in %s.',files(i).name);
    end
    if numel(s.window_position) < 2
        warning('SleepScan:singleWindow','%s has only one window; excluded from across-fly progression.',s.subject);
        continue
    end
    sleep(:,:,i) = interp1(s.window_position,s.sleep_mean.',grid,'linear').';
    randomSample(:,i) = get_random_sample_mean(s);
    included(i) = true;
end
if ~any(included)
    warning('SleepScan:noProgression','No flies have more than one window. Per-fly figures are still available.');
    return
end
group.subjects = cellfun(@(s) s.subject,summaries(included),'UniformOutput',false);
group.window_position = grid;
group.channel_list = channels;
group.sleep_by_fly = sleep(:,:,included);
group.random_sample_by_fly = randomSample(:,included);
group.nFlies = sum(isfinite(group.sleep_by_fly),3);
group.mean_sleep = mean(group.sleep_by_fly,3,'omitnan');
group.mean_random_sample = mean(group.random_sample_by_fly,2,'omitnan');
group.sem_sleep = std(group.sleep_by_fly,0,3,'omitnan')./sqrt(group.nFlies);
group.sem_random_sample = std(group.random_sample_by_fly,0,2,'omitnan')./sqrt(sum(isfinite(group.random_sample_by_fly),2));
group.sem_sleep(group.nFlies < 2) = NaN;
group.sem_random_sample(sum(isfinite(group.random_sample_by_fly),2) < 2) = NaN;
group.description = 'Equal fly weighting; linear interpolation from first to last possible window, not elapsed time. SEM across flies, not overlapping windows.';
stem = fullfile(outputDir,['All_flies_',first.condition,'_all_electrodes']);
activityTag = get_activity_tag(first);
heading = sprintf('%s: mean +/- SEM across %d flies',activityTag,nnz(included));
if nnz(included) == 1
    heading = sprintf('%s: 1 available fly (SEM unavailable)',activityTag);
end
draw_panels(100*grid,group.mean_sleep,group.sem_sleep,group.mean_random_sample,channels,heading, ...
    'Relative window position (0 = first, 100 = last; not elapsed time)',stem);
save([stem,'.mat'],'group','-v7.3');
end

function draw_panels(x,sleep,semSleep,randomSampleMean,channels,heading,xlabelText,stem)
fig = figure('Color','w','Position',[100 100 1500 850],'Visible','off');
cleanup = onCleanup(@() close(fig));
t = tiledlayout(fig,ceil(numel(channels)/5),5,'TileSpacing','compact','Padding','compact');
sleepColour = [0.85 0.4 0.12];
randomColour = [0.15 0.15 0.15];
for e = 1:numel(channels)
    ax = nexttile(t); hold(ax,'on');
    if ~isempty(semSleep)
        shade(ax,x,sleep(e,:),semSleep(e,:),sleepColour);
    end
    h1 = plot(ax,x,sleep(e,:),'Color',sleepColour,'LineWidth',0.4,'Marker','.','MarkerSize',4);
    if ~isempty(randomSampleMean) && isfinite(randomSampleMean(e))
        h2 = yline(ax,randomSampleMean(e),':','Color',randomColour,'LineWidth',0.6);
    else
        h2 = gobjects(1);
    end
    yline(ax,0,':','Color',[0.4 0.4 0.4],'LineWidth',0.3);
    if numel(x)>1, xlim(ax,[x(1),x(end)]); end
    title(ax,sprintf('E%d',channels(e)));
    set(ax,'Box','off','TickDir','out','LineWidth',0.3,'FontName','Arial', ...
        'FontSize',9,'XGrid','off','YGrid','off');
    if e == 1
        if isgraphics(h2)
            lg = legend(ax,[h1 h2],{'Sliding-window MI','Random-sample MI'}, ...
                'Orientation','horizontal','Box','off');
        else
            lg = legend(ax,h1,{'Sliding-window MI'}, ...
                'Orientation','horizontal','Box','off');
        end
        lg.Layout.Tile = 'north';
    end
end
title(t,heading,'Interpreter','none');
xlabel(t,xlabelText); ylabel(t,'Mean across ERP time points (bits)');
savefig(fig,[stem,'.fig']);
print(fig,[stem,'.png'],'-dpng','-r200');
print(fig,[stem,'.svg'],'-dsvg','-painters');
fprintf('Saved %s (.fig/.png/.svg)\n',stem);
end

function shade(ax,x,m,se,colour)
if all(~isfinite(se)), return; end
fill(ax,[x fliplr(x)],[m-se fliplr(m+se)],colour, ...
    'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
end

function activityTag = get_activity_tag(summary)
if isfield(summary,'activity_tag')
    activityTag = summary.activity_tag;
elseif isfield(summary,'activity_tag_2')
    activityTag = summary.activity_tag_2;
else
    activityTag = 'selected activity';
end
end

function randomSampleMean = get_random_sample_mean(summary)
if isfield(summary,'random_sample_mean')
    randomSampleMean = summary.random_sample_mean(:);
else
    randomSampleMean = nan(numel(summary.channel_list),1);
end
end
