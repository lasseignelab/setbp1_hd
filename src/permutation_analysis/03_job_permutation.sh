#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=96GB
#SBATCH --time=10:00:00
#SBATCH --partition=short

#### DATASETS ####

cap_container -c singularity "lizzyr/lw_dea:0.3.1"

singularity exec --cleanenv \
                 --containall \
                 -B "${LASSEIGNE_LAB_PATH}" \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/lw_dea_0.3.1.sif \
                  Rscript --vanilla "${CAP_PROJECT_PATH}"/src/permutation_analysis/03_postnatal_development.R \
                  -i "${CAP_DATA_PATH}/" \
                  -o "${CAP_RESULTS_PATH}"/permutation_analysis/
