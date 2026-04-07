#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=short
#SBATCH --time=00:30:00

#### DATASETS ####

cap_container -c singularity "lizzyr/lw_dea:0.3.1"
script="${CAP_PROJECT_PATH}"/src/permutation_analysis/04_prioritize_targets.R

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/lw_dea_0.3.1.sif \
                  Rscript --vanilla "${script}" \
                  -i "${CAP_DATA_PATH}/" \
                  -o "${CAP_RESULTS_PATH}"/permutation_analysis/
