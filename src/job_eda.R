#!/bin/bash
## run the 01_exploratory_data_analysis.R script and schedule this job with
## `sbatch 01_job-preprocessed_EDA.sh`
##

#SBATCH --job-name EDA-postproc   ## name that will show up in the queue
#SBATCH --ntasks=1  ## number of tasks (analyses) to run, i.e. # of nodes to request resource from
#SBATCH --cpus-per-task=12  ## the number of threads allocated to each task
#SBATCH --mem=64G  # memory per CPU core
#SBATCH --partition=short  ## the partition to run in (short == 12h max run time)
#
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-4

## Load modules
ml Singularity

## Variables
wd="$USER_DATA/230815_LW_ExerciseAD"

## Input list of datasets
dataset_list=(`cat ${wd}/data/dataset_names.txt`)
dataset=`echo ${dataset_list[$SLURM_ARRAY_TASK_ID]}`

## Code
cd ${wd}
singularity exec -B ${wd} ${wd}/bin/docker/lw_eda_1.0.0.sif Rscript --vanilla ${wd}/src/01_exploratory_data_analysis.R \
-r ${wd}/data/processed_rse/${dataset}_rse.rds \
-t "sra_attribute.exercise" \
-o ${wd}/results/eda/processed/
