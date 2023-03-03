# CoI-pipeline
Pipeline for estimating mutual information and co_information using gaussian copula estimation

This pipeline uses GCMI (Ince et al, 2016) to estimate mutual information and co-information from electrophysiological data of auditory oddball paradigms. 

You will also need GCMI-functions and CMI-functions for the pipeline to work.

Data folder shows the folder structure the data has to be organised in for the pipeline to work. Empty text-files show the naming conventions to use.

SCRIPTS FOR THE PIPELINE 

Get_param.m:

Defines all of the attributes of the data, and has to be changed according to the analysed data. Only change the parameters here to run the pipeline with your
data

main.m:

Function build to function with the Cambridge HPC, takes in a SLURM-job as an input.

ImpiEEG.m:

Imports and preprocesses data-based on parameters defined in Get_param.m

EoI_data.m:

Imports electrodes of interests EoI.mat-files from the Data/EoI -folder

Get_COI.m:

Estimates Co-information for the specified electrode pair (decided by the slurm-job id main takes as input)

MI_ERP:

Estimates mutual information and performs parametric permutation testing. Saves significant electrodes into a mat-file. (Run locally and upload said mat-file into DATA/EOI
