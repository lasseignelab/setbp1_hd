#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=short

cap_container -c singularity "projectassistant/sitemap-scraper:latest"

singularity exec --cleanenv \
                 --containall \
                 -B "${CAP_PROJECT_PATH}" \
                 "${CAP_CONTAINER_PATH}"/sitemap-scraper_latest.sif \
                  python3 "${CAP_PROJECT_PATH}"/src/13_drug_safety.py \
                  -i "${CAP_RESULTS_PATH}"/signature_reversion/brain_drugs_filtered.csv \
                  -o "${CAP_RESULTS_PATH}"/signature_reversion
