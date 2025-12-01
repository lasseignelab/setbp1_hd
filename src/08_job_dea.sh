#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=express
#SBATCH --array=0-3

#### DATASETS ####
dataset_names=("cardo" "shaw" "tanaka1" "tanaka2")
dataset=${dataset_names[$SLURM_ARRAY_TASK_ID]}

cap_container -c singularity "lizzy/lw_dea:0.1.1"

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/lw_dea_0.1.1.sif \
                  Rscript --vanilla "${CAP_PROJECT_PATH}"/src/08_find_degs.R \
                  -r "${CAP_DATA_PATH}"/processed/"${dataset}"_prepped_rse.rds \
                  -s "${dataset}" \
                  -o "${CAP_PROJECT_PATH}"/results/dea/"${dataset}"
