#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=express

cap_container -c singularity "lizzyr/sigsearch:1.0.0"

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 -B "${LASSEIGNE_LAB_PATH}" \
                 "${CAP_CONTAINER_PATH}"/sigsearch_1.0.0.sif \
                  Rscript --vanilla "${CAP_PROJECT_PATH}"/src/figures/figure_5.R \
                  -i "${CAP_RESULTS_PATH}"/signature_reversion/ \
                  -d "${CAP_DATA_PATH}"/dailymed/ \
                  -o "${CAP_RESULTS_PATH}"/figures/
