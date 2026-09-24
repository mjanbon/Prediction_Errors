# Sleep-trial window scan

The original E14 figure in `Data/MI_Data/tmp_R060721/sleep_window_scan`
comes from `run_R060721_sleep_window_MI_comparison.m`. That helper is unchanged.
`run_sleep_window_MI_comparison.m` extends that analysis to all 15 retained
electrodes and all participants using the single activity selected in
`Max_get_param.m`.

## Run locally

Add `Helpers` to the MATLAB path, then run:

```matlab
run_sleep_window_MI_comparison()       % all configured flies
run_sleep_window_MI_comparison(0,1)    % first fly only
run_sleep_window_MI_comparison(0,[],struct('windowFraction',0.5)) % sliding range in every fly
run_sleep_window_MI_comparison(0,[],struct('windowFraction',0.1,'windowStep',50)) % local 10% windows
```

The default is a descriptive scan with no permutations. This does not change
the MI estimates used for these figures. Optional settings can be supplied:

```matlab
options = struct('channels',15:-1:1,'windowTrials',[], ...
                 'windowStep',[],'nPerm',0);
run_sleep_window_MI_comparison(0,[],options)
```

An empty window length uses the smaller available deviant/standard count for
the activity selected by `Max_get_param`, multiplied by `windowFraction`
(default 1). Use 0.5 for a denser progression in flies with enough selected
activity trials. These runs are saved separately in a `windowFraction_0p5`
subfolder. An
explicit `windowTrials` overrides the fraction and uses a `trials_N` subfolder.
The default step is one tenth of the chosen length. Set `windowStep` directly
for more windows; for example `windowStep=50` advances the window by 50 trials
at a time and saves into a `step_50` subfolder. The current Apocrita script
uses `windowFraction=0.1`, so each window contains 10% of that fly's balanced
deviant/standard trial count. Each window contains N deviant and N standard
trials from the selected activity. `nPerm` is currently ignored because this
scan no longer calculates a wake-minus-sleep comparison.

Outputs go under:

```text
Data/MI_Data/sleep_window_scan_all_electrodes/<activity_tag>_BSLEEP/
```

Each fly has a numerical MAT summary and a 15-panel FIG/PNG/SVG figure.
Orange is MI for the selected activity, averaged over the full 100-sample ERP.
The across-fly figure gives each fly equal weight
and shows SEM across flies. It interpolates each curve to 101 positions from
the first to last possible window. A fly with only one possible window has a
per-fly figure but is excluded from the progression average. The saved group
MAT records included subjects, each fly's interpolated curve and sample counts.

Replot or aggregate existing summaries without loading LFP data:

```matlab
plot_sleep_window_MI_scan(outputDir, 'R060721') % one fly
plot_sleep_window_MI_scan(outputDir)            % across flies
```

## Run on Apocrita

Upload the three new `.m` helpers and submit
`Code/apocrita_slurm_submit_sleep_window_MI_all_electrodes.sh`.
Its array tasks 1-8 select the configured participants, one per task, using
`windowFraction=0.5` to allow sliding in every fly. After
all tasks finish, run the aggregation command at the bottom of that script,
or copy the summary MAT files locally and call `plot_sleep_window_MI_scan`.
The script uses the project's existing log directory and hardcoded cluster paths.

## Interpretation and compatibility

- This scans **colour-concatenated trial indices**, not elapsed time, sleep
  bouts, or time within a bout. The existing loader's colour ordering is
  retained for comparability with the E14 figure. Although the complete
  loaded condition is colour-balanced, an individual window need not be.
- Windows overlap; they are not independent samples for statistical testing.
  The across-fly SEM uses flies, not windows, as its observations.
- Electrode labels use the existing E1-E15 convention after removing the
  reference channel. Original channel numbers are saved. A common valid
  channel mask is used for deviant and standard data, and recordings with
  further channel losses are skipped instead of silently shifting the
  electrode map.
- Local files with 20 rather than 100 samples are skipped. Use the correctly
  sampled LFP files on the cluster for those recordings. Missing flies are
  reported and are not silently represented in the group average.
- The separate `sleep_window_length_scan` analysis varies window length;
  this implementation varies window position at a fixed length within fly.
