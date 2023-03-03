Pipeline for estimating MI & CoI with GCMI
=====================================================
This repository contains scripts for estimating MI and CoI using gaussian copula estimation ([Ince, 2016](https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.23471)). Designed to be run on the Cambridge HPC-server.

I have tested this with the marmoset-data I have, and it is working for me at least, but let me know of any problems!

Data-folder
----------------------------

The data-folder contains the folder structure you should use with pipeline for it work. 

Empty text files show the naming convention for the data-files: 'Participant_Condition_sta.set'

EoI-folder should contain the mat-files with the electrodes of interest. 


Get MI and electrodes of interest
----------------------------
```
MI_ERP.m
```
RUN THIS FIRST LOCALLY. Estimates MI and performs parametric permutation testing. Saves channels with significant MI to a mat-file (that you will have to upload to the EoI-folder for the pipeline to work).

HPC Pipeline
----------------------------
```
Get_param.m
```
If data is saved in the right folder structure/naming conventions (and the code functions as I've intended) this is the only function you will have to edit. Change all the parameters according to how you want to analyse your data.

```
main.m
```
Should function with the Cambridge HPC. Takes in a slurm-job as an input, which specifies the electrode pair MI/CoI will be estimated for.

```
EoI_data.m
```
Loads the electrodes of interest from the mat-files in the EoI-folder

```
ImpiEEG.m
```
Loads the data onto a matlab-structure based on the parameters specified in Get_param.

```
Get_COI.m
```
Estimates CoI for the specified electrode pair (which is decided by the slurm-job index that the main.m takes as input)
