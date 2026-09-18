#!/bin/bash
#SBATCH -o /data/SBBS-PIDProject/Maxime/Logs/Drosophila_MI/MI_LFP_peak_summary_%j.log
#SBATCH --mem-per-cpu=16G
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH --partition=compute

module load matlab

CODE_DIR="/data/SBBS-PIDProject/Maxime/Prediction_Errors/CoI-pipeline/Code"
PLOT_DIR="/data/SBBS-PIDProject/Maxime/Prediction_Errors/CoI-pipeline/Plotting/Refactored_Plotting"

matlab -nodisplay -nosplash -r "addpath('${CODE_DIR}'); addpath('${PLOT_DIR}'); plot_MI_LFP_peak_summary(2); quit"
