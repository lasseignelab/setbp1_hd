#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=express

cap_container -c singularity "lizzyr/lw_dea:0.3.1"

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/lw_dea_0.3.1.sif \
                  Rscript --vanilla "${CAP_PROJECT_PATH}"/src/figures/figure_3.R \
                  -i "${CAP_RESULTS_PATH}"/consensus_signature/ \
                  -d "${CAP_RESULTS_PATH}"/dea/ \
                  -o "${CAP_RESULTS_PATH}"/figures/
