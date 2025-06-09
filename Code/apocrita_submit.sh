#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -j y
#$ -l h_vmem=8G
#$ -o /data/SBBS-PIDProject/Maxime/Logs/CoI/
#$ -e /data/SBBS-PIDProject/Maxime/Logs/CoI/
#$ -l h_rt=239:00:00
#$ -t 1-1800


# Load MATLAB ready for processing.
module load matlab


matlab -nodisplay -nosplash -r "max_main(${SGE_TASK_ID}); quit"