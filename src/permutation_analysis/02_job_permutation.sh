#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --partition=short

#### DATASETS ####

cap_container -c singularity "lizzyr/lw_dea:0.3.1"

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/lw_dea_0.3.1.sif \
                  Rscript --vanilla "${CAP_PROJECT_PATH}"/src/permutation_analysis/02_spatiotemporal_perm.R \
                  -i "${CAP_DATA_PATH}/" \
                  -o "${CAP_RESULTS_PATH}"/permutation_analysis/
