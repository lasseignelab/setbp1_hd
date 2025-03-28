#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=express
#SBATCH --array=0-2


#### FASTQ DIRECTORIES ####
dataset_names=("cardo" "shaw" "tanaka")
dataset=${dataset_names[$SLURM_ARRAY_TASK_ID]}
fastq_dir="$CAP_DATA_PATH/raw/$dataset"

#### CREATE SAMPLESHEETS ####
script="$CAP_PROJECT_PATH/src/03_fastq_dir_to_samplesheet.py"
outfile="$CAP_DATA_PATH/samplesheet_${dataset}.csv"

# create samplesheet
python3 "$script" "$fastq_dir" "$outfile" \
    -r1 "_1.fastq.gz" \
    -r2 "_2.fastq.gz"

echo "Samplesheet saved to $outfile."
