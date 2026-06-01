#!/bin/bash
#SBATCH -o /data/SBBS-PIDProject/Maxime/Logs/Window_Length_Analysis/sleep_window_length_scan_%j.log
#SBATCH --mem-per-cpu=16G    # Request 16G per cpu
#SBATCH -t 239:00:00         # Request 239 hours runtime
#SBATCH -n 1                # Request 1 core
#SBATCH -a 1-8           # Array Job-ID 1-8 (adjust to number of subjects)
#SBATCH --partition=compute # Request the compute partition

# Load MATLAB module
module load matlab

matlab -nodisplay -nosplash -r "run_sleep_window_length_scan_array(${SLURM_ARRAY_TASK_ID}); quit"
