#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=medium


#### ENVIRONMENT ####
# cheaha profile info https://github.com/nf-core/configs/blob/master/docs/cheaha.md
ml Anaconda3
conda activate nfcore-env

# Run pipeline
nextflow run nf-core/rnaseq -r 3.18.0 \
    --input "$CAP_DATA_PATH/_samplesheet.csv" \
    --outdir "$CAP_DATA_PATH" \
    -profile cheaha \
    --aligner star_salmon \
    --pseudo_aligner salmon \
    --fasta "$CAP_DATA_PATH/GRCh38.primary_assembly.genome.fa" \
    --gtf "$CAP_DATA_PATH/gencode.v47.primary_assembly.annotation.gtf" \
    --gencode TRUE \
    --save_reference TRUE
