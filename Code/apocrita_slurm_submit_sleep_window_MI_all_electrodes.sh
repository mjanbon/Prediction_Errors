#!/bin/bash
#SBATCH -o /data/SBBS-PIDProject/Maxime/Logs/Window_Length_Analysis/sleep_window_MI_%A_%a.log
#SBATCH --mem=32G
#SBATCH -t 239:00:00
#SBATCH -n 1
#SBATCH -a 1-8
#SBATCH --partition=compute

module load matlab
HELPER_DIR="/data/SBBS-PIDProject/Maxime/Prediction_Errors/CoI-pipeline/Helpers"

# One fly per task. Half the matched length allows sliding even in flies
# with fewer sleep than wake trials. Set windowFraction=1 for the old length.
matlab -batch "addpath('${HELPER_DIR}'); run_sleep_window_MI_comparison(2,${SLURM_ARRAY_TASK_ID},struct('windowFraction',0.5));"

# After all eight tasks finish, aggregate their saved summaries in MATLAB:
# plot_sleep_window_MI_scan('/data/SBBS-PIDProject/Maxime/Drosophila_Data/MI_Data/sleep_window_scan_all_electrodes/sleep_BSLEEP/windowFraction_0p5')
