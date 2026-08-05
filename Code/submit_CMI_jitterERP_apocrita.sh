#!/bin/bash
#SBATCH -o /data/SBBS-PIDProject/Maxime/Logs/Jitter/CMI_jitterERP_%A_%a.log
#SBATCH --mem-per-cpu=32G
#SBATCH -t 239:00:00
#SBATCH -n 1
#SBATCH -a 1-7
#SBATCH --partition=compute

# Submit from Apocrita with:
#   sbatch submit_CMI_jitterERP_apocrita.sh
#
# Each array task analyses one jitterERP entry:
#   SLURM_ARRAY_TASK_ID=1 -> jitterERP(1)
#   ...
#   SLURM_ARRAY_TASK_ID=7 -> jitterERP(7)

module load matlab

CODE_DIR="/data/SBBS-PIDProject/Maxime/Prediction_Errors/CoI-pipeline/Code"
JITTER_ERP_FILE="/data/SBBS-PIDProject/Maxime/Drosophila_Data/Drosophila_LFP/Rdaytime_jitteringExperiment_jitteringExperiment_jitterERP.mat"
OUTPUT_DIR="/data/SBBS-PIDProject/Maxime/Drosophila_Data/Drosophila_LFP/CMI_Data"
N_PERM=2

matlab -nodisplay -nosplash -r "addpath('${CODE_DIR}'); main_CMI_jitterERP('${JITTER_ERP_FILE}','${OUTPUT_DIR}',${SLURM_ARRAY_TASK_ID},${N_PERM}); quit"
