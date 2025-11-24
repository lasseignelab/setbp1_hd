#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=medium

cap_container -c singularity "lizzyr/sigsearch:0.0.0"

echo "Keeps going when the container exists."
