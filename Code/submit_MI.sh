#!/bin/bash
#SBATCH -p icelake
#SBATCH -t 00:05:00
#SBATCH -N = 1
#SBATCH -J Drosophila_MI_R28
#SBATCH -o logs/Drosophila_MI_R28_%A_%a.out
#SBATCH -e logs/Drosophila_MI_R28_%A_%a.err
#SBATCH --mem 24000

###############################
#  DO NOT CHANGE THESE LINES
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
###############################
cd /home/mj649/rds/hpc-work/
module load matlab
cd /home/mj649/rds/hpc-work/Prediction_Errors/Code
matlab -nodesktop -nosplash -r "main_MI_ERP(1); quit"