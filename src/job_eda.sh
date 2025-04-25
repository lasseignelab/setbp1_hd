#!/bin/bash

#SBATCH --ntasks=1  ## number of tasks (analyses) to run, i.e. # of nodes to request resource from
#SBATCH --cpus-per-task=12  ## the number of threads allocated to each task
#SBATCH --mem=64G  # memory per CPU core
#SBATCH --partition=short  ## the partition to run in (short == 12h max run time)
#SBATCH --array=0-1

#### DATASETS ####
dataset_names=("cardo" "shaw" "tanaka")
dataset=${dataset_names[$SLURM_ARRAY_TASK_ID]}

#### INPUTS ####
counts="${CAP_DATA_PATH}/${dataset}/star_salmon/salmon.merged.gene_counts.SummarizedExperiment.rds"
metadata="${CAP_PROJECT_PATH}/doc/${dataset}_metadata.csv"


## Input list of datasets
#dataset_list=(`cat ${wd}/data/dataset_names.txt`)
#dataset=`echo ${dataset_list[$SLURM_ARRAY_TASK_ID]}`

#### RUN EDA SCRIPT ####
singularity exec -B ${CAP_PROJECT_PATH} ${CAP_CONTAINER_PATH}/lw_eda_2.1.1.sif Rscript --vanilla ${CAP_PROJECT_PATH}/src/eda.R \
-r ${counts} \
-m ${metadata} \
-o ${CAP_PROJECT_PATH}/results/eda/
