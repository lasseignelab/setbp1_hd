#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=medium


#### ENVIRONMENT ####
# cheaha profile info https://github.com/nf-core/configs/blob/master/docs/cheaha.md
ml Anaconda3
conda activate nfcore-env

# nextflow memory requirements
NXF_OPTS='-Xms1g -Xmx4g'echo "Nextflow environment settings:"
env | grep "NXF"

# Run pipeline
nextflow run nf-core/rnaseq -r 3.18.0 \
    --input "$CAP_DATA_PATH/raw/samplesheet/all_samplesheet.csv" \
    --outdir "$CAP_DATA_PATH" \
    -profile cheaha \
    --aligner star_salmon \
    --pseudo_aligner salmon \
    --fasta "$CAP_DATA_PATH/GRCh38.primary_assembly.genome.fa" \
    --gtf "$CAP_DATA_PATH/gencode.v45.annotation.gtf" \
    --gencode TRUE \
    --save_reference TRUE
