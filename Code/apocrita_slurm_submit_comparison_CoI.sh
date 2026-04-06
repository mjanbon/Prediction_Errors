#!/bin/bash
#SBATCH -o /data/SBBS-PIDProject/Maxime/Logs/CoI_Comparison/CoI_permute_comparison_%j.log
#SBATCH --mem-per-cpu=16G    # Request 16G per cpu
#SBATCH -t 239:00:00         # Request 239 hours runtime
#SBATCH -n 1                # Request 1 core
#SBATCH -a 3-1800               # Array Job-ID 3-1800
#SBATCH --partition=compute # Request the compute partition


# Load MATLAB ready for processing.
module load matlab


matlab -nodisplay -nosplash -r "max_main_comparison(${SLURM_ARRAY_TASK_ID}); quit"