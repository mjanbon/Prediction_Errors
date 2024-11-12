#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -j y
#$ -l h_vmem=4G
#$ -o ~
#$ -l h_rt=1:0:0
#$ -t 1-225


# Load MATLAB ready for processing.
module load matlab


matlab -nodisplay -nosplash -r "max_main(${SGE_TASK_ID}); quit"