#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -j y
#$ -l h_vmem=64G
#$ -o /data/SBBS-PIDProject/Maxime/Logs/Drosophila_MI/
#$ -e /data/SBBS-PIDProject/Maxime/Logs/Drosophila_MI/
#$ -l h_rt=239:0:0
#$ -l rocky


# Load MATLAB ready for processing.
module load matlab

matlab -nodisplay -nosplash -r "main_MI_ERP_comparison(1); quit"