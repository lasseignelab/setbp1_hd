#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=medium


#### VARIABLES ####
ENV_NAME="nfcore-env"
dataset="tanaka"
fastq_dir="$CAP_DATA_PATH/raw/$dataset"
out_dir="$CAP_DATA_PATH/$dataset"
run_dir="$USER_SCRATCH/nfcore_runs/${dataset}"
work_dir="$USER_SCRATCH/nfcore_work/${dataset}"

mkdir -p "$run_dir" "$out_dir" "$work_dir"
cd "$run_dir"

#### ENVIRONMENT ####
ml Anaconda3/2023.07-2
conda activate $ENV_NAME

# Run pipeline --cheaha profile info https://github.com/nf-core/configs/blob/master/docs/cheaha.md
nextflow run nf-core/rnaseq -r 3.18.0 \
    --input "$CAP_DATA_PATH/samplesheet_${dataset}.csv" \
    --outdir "$out_dir" \
    -profile cheaha \
    --aligner star_salmon \
    --pseudo_aligner salmon \
    --fasta "$CAP_DATA_PATH/GRCm39.primary_assembly.genome.fa" \
    --gtf "$CAP_DATA_PATH/gencode.vM31.primary_assembly.annotation.gtf" \
    --gencode TRUE \
    --save_reference TRUE
