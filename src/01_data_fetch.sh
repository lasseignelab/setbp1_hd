#!/usr/bin/env bash

# SLURM settings
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=medium

#### ENVIRONMENT ####
# cheaha profile info https://github.com/nf-core/configs/blob/master/docs/cheaha.md
ml Anaconda3
conda activate nfcore-env

# nextflow memory requirements
NXF_OPTS='-Xms1g -Xmx4g'

# directory to output data
out_dir="$CAP_PROJECT_PATH/data/raw"
mkdir -p "$out_dir"

nextflow run nf-core/fetchngs \
   -profile cheaha \
   --download_method sratools \
   --input "$CAP_PROJECT_PATH/doc/dataset_ids.csv" \
   --nf_core_pipeline rnaseq \
   --outdir "$out_dir"
