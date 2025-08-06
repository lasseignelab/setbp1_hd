#!/usr/bin/env bash

#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=12G
#SBATCH --time=72:00:00
#SBATCH --partition=long
#SBATCH --array=0-1

#### VARIABLES ####
ENV_NAME="nfcore-env"
dataset_names=("cardo" "shaw")
dataset="${dataset_names[$SLURM_ARRAY_TASK_ID]}"
out_dir="$CAP_DATA_PATH/$dataset"
run_dir="$USER_SCRATCH/nfcore_runs/${dataset}"
work_dir="$USER_SCRATCH/nfcore_work/${dataset}"

# create isolated run and output directories to run in parallel
mkdir -p "$run_dir" "$out_dir" "$work_dir"
cd "$run_dir" || { echo "Error: Unable to CD to ${run_dir}" >&2; exit 1; }

#### ENVIRONMENT ####
ml Anaconda3/2023.07-2
conda activate $ENV_NAME

# run pipeline for each data set
nextflow run nf-core/rnaseq -r 3.18.0 \
    --input "$CAP_DATA_PATH/raw/samplesheet_${dataset}.csv" \
    --outdir "$out_dir" \
    -profile cheaha \
    -work-dir "$work_dir" \
    --aligner star_salmon \
    --pseudo_aligner salmon \
    --fasta "$CAP_DATA_PATH/GRCh38.primary_assembly.genome.fa" \
    --gtf "$CAP_DATA_PATH/gencode.v47.primary_assembly.annotation.gtf" \
    --gencode TRUE
