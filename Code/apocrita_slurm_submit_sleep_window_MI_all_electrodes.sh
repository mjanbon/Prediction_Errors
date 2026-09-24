#!/bin/bash
#SBATCH -o /data/SBBS-PIDProject/Maxime/Logs/Window_Length_Analysis/sleep_window_MI_%A_%a.log
#SBATCH --mem=32G
#SBATCH -t 239:00:00
#SBATCH -n 1
#SBATCH -a 1-8
#SBATCH --partition=compute

module load matlab
HELPER_DIR="/data/SBBS-PIDProject/Maxime/Prediction_Errors/Helpers"

# One fly per task. Each window uses 10% of that fly's balanced trial count.
# windowStep=50 advances the window by 50 trials for a dense local scan.
matlab -batch "addpath('${HELPER_DIR}'); run_sleep_window_MI_comparison(2,${SLURM_ARRAY_TASK_ID},struct('windowFraction',0.1,'windowStep',50));"

# After all eight tasks finish, aggregate their saved summaries in MATLAB:
# plot_sleep_window_MI_scan('/data/SBBS-PIDProject/Maxime/Drosophila_Data/MI_Data/sleep_window_scan_all_electrodes/sleep_BSLEEP/windowFraction_0p1/step_50')
