#!/bin/bash

#SBATCH --ntasks=1  ## number of tasks (analyses) to run, i.e. # of nodes to request resource from
#SBATCH --cpus-per-task=4  ## the number of threads allocated to each task
#SBATCH --mem=12G  # memory per CPU core
#SBATCH --partition=express  ## the partition to run in (short == 12h max run time)
#SBATCH --array=0-2

#### DATASETS ####
dataset_names=("cardo" "shaw" "tanaka")
dataset=${dataset_names[$SLURM_ARRAY_TASK_ID]}

#### INPUTS ####
counts="${CAP_DATA_PATH}/${dataset}/star_salmon/salmon.merged.gene_counts.SummarizedExperiment.rds"
metadata="${CAP_PROJECT_PATH}/doc/${dataset}_metadata.csv"

#### RUN EDA SCRIPT ####
singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 -B /data/project/lasseigne_lab \
                 "${CAP_CONTAINER_PATH}/lw_eda_2.1.1.sif" Rscript --vanilla "${CAP_PROJECT_PATH}/src/05_data_cleaning.R" \
                 -f "${CAP_PROJECT_PATH}/src/functions.R" \
                 -c "${counts}" \
                 -m "${metadata}" \
                 -n "${dataset}" \
                 -o "${CAP_DATA_PATH}/processed/"
