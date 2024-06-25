#!/bin/bash
#SBATCH -p icelake
#SBATCH -t 00:10:00
#SBATCH --array=225-450
#SBATCH -J Drosophila_CoI_R23
#SBATCH -o logs/Drosophila_CoI_R23_%A_%a.out
#SBATCH -e logs/Drosophila_CoI_R23_%A_%a.err
#SBATCH --mem 24000

###############################
#  DO NOT CHANGE THESE LINES
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
###############################

module load matlab
cd /home/mj649/rds/hpc-work/Prediction_Errors/Code
matlab -nodesktop -nosplash -r "max_main(${SLURM_ARRAY_TASK_ID}); quit"
