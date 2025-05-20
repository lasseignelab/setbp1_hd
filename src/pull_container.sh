#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=medium

cap_container -c singularity "lizzyr/lw_dea:0.1.1"

echo "Keeps going when the container exists."
