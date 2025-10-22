#!/bin/bash

#SBATCH --ntasks=1  ## number of tasks (analyses) to run, i.e. # of nodes to request resource from
#SBATCH --cpus-per-task=12  ## the number of threads allocated to each task
#SBATCH --mem=64G  # memory per CPU core
#SBATCH --partition=express ## the partition to run in (short == 12h max run time)

cap_container -c singularity "lizzyr/lw_models:0.1.1"

#### RUN EDA SCRIPT ####
singularity exec --cleanenv \
                 --containall \
                 -B /data/project/lasseigne_lab \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}/lw_models_0.1.1.sif" Rscript --vanilla "${CAP_PROJECT_PATH}/src/07_combined_eda.R" \
                 -f "${CAP_PROJECT_PATH}/src/functions.R" \
                 -r "${CAP_DATA_PATH}/processed/" \
                 -o "${CAP_PROJECT_PATH}/results/eda/combined_datasets/"
