#!/bin/bash
#SBATCH -o /data/SBBS-PIDProject/Maxime/Logs/SLURM_Test/test_slurm_%j.log
#SBATCH -e /data/SBBS-PIDProject/Maxime/Logs/SLURM_Test/test_slurm_%j.err
#SBATCH --mem-per-cpu=2G    # Request 2G per cpu
#SBATCH -t 00:05:00         # Request 5 mins runtime
#SBATCH -n 1                # Request 1 core

# Load MATLAB module
module load matlab

# Run the test script
matlab -nodisplay -nosplash -r "test_slurm; quit"
