#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=express

#### CONDA ENVIRONMENT ####
YML_FILE="$CAP_CONDA_PATH/sra-env.yml"
ENV_NAME="sra-env"

ml Anaconda3/2023.07-2

# Check if the environment already exists
if conda info --envs | grep -q "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists."
else
    echo "Creating environment '$ENV_NAME' from $YML_FILE..."
    conda env create -f "$YML_FILE"
fi
