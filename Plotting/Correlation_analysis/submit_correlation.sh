#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -j y
#$ -l h_vmem=8G
#$ -o /data/SBBS-PIDProject/Maxime/Logs/Correlation/
#$ -e /data/SBBS-PIDProject/Maxime/Logs/Correlation/
#$ -l h_rt=1:00:00
#$ -t 3-840

# Load MATLAB ready for processing.
module load matlab

# Run the correlation analysis with task ID
matlab -nodisplay -nosplash -r "correlation_main(${SGE_TASK_ID}); quit"