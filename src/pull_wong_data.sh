#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=medium
#SBATCH --partition=express

cap_data_download \
  --subdirectory "wong" \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-64074-x/MediaObjects/41467_2025_64074_MOESM14_ESM.zip"

unzip -d "${CAP_DATA_PATH}/wong" "${CAP_DATA_PATH}/wong/41467_2025_64074_MOESM14_ESM.zip"
