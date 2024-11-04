#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -j y
#$ -l h_vmem=4G
#$ -o ~
#$ -l h_rt=1:0:0


# Load MATLAB ready for processing.
module load matlab

matlab -nodisplay -nosplash -r "main_MI_ERP(1); quit"