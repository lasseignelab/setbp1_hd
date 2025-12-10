#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=express
#SBATCH --array=0-1

#### DATASETS ####
dataset_names=("cardo" "shaw")
dataset=${dataset_names[$SLURM_ARRAY_TASK_ID]}

cap_container -c singularity "lizzyr/lw_dea:0.3.1"

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/lw_dea_0.3.1.sif \
                  Rscript --vanilla "${CAP_PROJECT_PATH}"/src/09_analyze_degs.R \
                  -d "${CAP_RESULTS_PATH}"/dea/"${dataset}"/significant_degs.csv \
                  -o "${CAP_RESULTS_PATH}"/dea/"${dataset}"/
