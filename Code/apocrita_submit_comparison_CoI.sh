#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -j y
#$ -l h_vmem=8G
#$ -o /data/SBBS-PIDProject/Maxime/Logs/CoI_Comparison/
#$ -e /data/SBBS-PIDProject/Maxime/Logs/CoI_Comparison/
#$ -l h_rt=239:00:00
#$ -t 1-2


# Load MATLAB ready for processing.
module load matlab


matlab -nodisplay -nosplash -r "max_main_comparison(${SGE_TASK_ID}); quit"