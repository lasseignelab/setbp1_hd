#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=medium


#### VARIABLES ####
YML_FILE="$CAP_PROJECT_PATH/bin/conda/nfcore-env.yml"
ENV_NAME="nfcore-env"
OUT_DIR="$CAP_DATA_PATH/tanaka/"

mkdir -p $OUT_DIR

#### ENVIRONMENT ####
# cheaha profile info https://github.com/nf-core/configs/blob/master/docs/cheaha.md
ml Anaconda3

# Check if the environment already exists
if conda info --envs | grep -q "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists."
else
    echo "Creating environment '$ENV_NAME' from $YML_FILE..."
    conda env create -f $YML_FILE
fi
conda activate $ENV_NAME

#### NF SETUP ####
# nextflow memory requirements
NXF_OPTS='-Xms1g -Xmx4g'echo "Nextflow environment settings:"
env | grep "NXF"

# Run pipeline
nextflow run nf-core/rnaseq -r 3.18.0 \
    --input "$CAP_DATA_PATH/samplesheet_tanaka.csv" \
    --outdir "$OUT_DIR" \
    -profile cheaha \
    --aligner star_salmon \
    --pseudo_aligner salmon \
    --fasta "$CAP_DATA_PATH/GRCm39.primary_assembly.genome.fa" \
    --gtf "$CAP_DATA_PATH/gencode.vM31.primary_assembly.annotation.gtf" \
    --gencode TRUE \
    --save_reference TRUE
