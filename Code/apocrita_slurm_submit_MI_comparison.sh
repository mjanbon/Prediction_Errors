#!/bin/bash
#SBATCH -o /data/SBBS-PIDProject/Maxime/Logs/Drosophila_MI/MI_permute_comparison_%j.log
#SBATCH --mem-per-cpu=64G    # Request 64G per cpu
#SBATCH -t 239:00:00         # Request 239 hours runtime
#SBATCH -n 1                # Request 1 core
#SBATCH --partition=compute # Request the compute partition

# Load MATLAB ready for processing.
module load matlab

matlab -nodisplay -nosplash -r "main_MI_ERP_comparison(1); quit"