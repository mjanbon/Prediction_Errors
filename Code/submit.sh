#!/bin/bash
#SBATCH -p icelake
#SBATCH -t 12:00:00
#SBATCH --array=1-144
#SBATCH -J VT_elecs
#SBATCH -o logs/VT_elec_%A_%a.out
#SBATCH -e logs/VT_elec_%A_%a.err
#SBATCH --mem 24000

###############################
#  DO NOT CHANGE THESE LINES
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
###############################

module load matlab
cd /home/jma201/coi-ieeg
matlab -nodesktop -nosplash -r "main_VT(${SLURM_ARRAY_TASK_ID}); quit"
