#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=express

cap_container -c singularity "lizzyr/sigsearch:0.0.0"

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 -B "${LASSEIGNE_LAB_PATH}" \
                 "${CAP_CONTAINER_PATH}"/sigsearch_0.0.0.sif \
                  Rscript --vanilla "${CAP_PROJECT_PATH}"/src/11_signature_reversion.R \
                  -i "${CAP_DATA_PATH}"/ \
                  -o "${CAP_RESULTS_PATH}"/
