#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=express

cap_container -c singularity "lizzy/lw_dea:0.3.1"

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/lw_dea_0.3.1.sif \
                  Rscript --vanilla "${CAP_PROJECT_PATH}"/src/10_consensus_signature.R \
                  -d "${CAP_RESULTS_PATH}"/dea/ \
                  -o "${CAP_RESULTS_PATH}"/consensus_signature/
