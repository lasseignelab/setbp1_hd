#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=long

#### VARIABLES ####
YML_FILE="$CAP_PROJECT_PATH/bin/conda/nfcore-env.yml"
ENV_NAME="nfcore-env"

dataset_names=("cardo" "shaw")

#### ENVIRONMENT ####
ml Anaconda3
conda activate $ENV_NAME

# Loop through each dataset to run pipeline sequentially 
#dataset=${dataset_names[$SLURM_ARRAY_TASK_ID]}
#fastq_dir="$CAP_DATA_PATH/raw/$dataset"
#out_dir="$CAP_DATA_PATH/$dataset"


for dataset in "${dataset_names[@]}"; do
    fastq_dir="$CAP_DATA_PATH/raw/$dataset"
    out_dir="$CAP_DATA_PATH/$dataset"

    mkdir -p $out_dir


    # Run pipeline -- cheaha profile info https://github.com/nf-core/configs/blob/master/docs/cheaha.md
    nextflow run nf-core/rnaseq -r 3.18.0 \
        --input "$CAP_DATA_PATH/samplesheet_${dataset}.csv" \
        --outdir "$out_dir" \
        -profile cheaha \
        --aligner star_salmon \
        --pseudo_aligner salmon \
        --fasta "$CAP_DATA_PATH/GRCh38.primary_assembly.genome.fa" \
        --gtf "$CAP_DATA_PATH/gencode.v47.primary_assembly.annotation.gtf" \
        --gencode TRUE
done
