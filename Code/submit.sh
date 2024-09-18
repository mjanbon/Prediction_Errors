#!/bin/bash
#SBATCH -p icelake
#SBATCH -t 00:30:00
#SBATCH --array=1-225
#SBATCH -J D_Start_Sleep_CoI_R06
#SBATCH -o logs/D_Start_Sleep_CoI_R06_%A_%a.out
#SBATCH -e logs/D_Start_Sleep_CoI_R06_%A_%a.err
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
matlab -nodesktop -nosplash -r "max_main(${SLURM_ARRAY_TASK_ID}); quit"
