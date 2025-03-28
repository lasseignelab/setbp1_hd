#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=short

# script for creating samples# script for creating samplesheet from directory of fastqs
file="$CAP_PROJECT_PATH/src/03_fastq_dir_to_samplesheet.py"
url="https://raw.githubusercontent.com/nf-core/rnaseq/refs/heads/master/bin/fastq_dir_to_samplesheet.py"

dataset=ADDHERE
fastq_dir="$CAP_DATA_PATH/raw/$dataset"

# check if file exists, pull if needed
if [! -f "$file"]; then
    echo "File $file not found. Downloading..."
    wget -P "$CAP_PROJECT_PATH/src" "$url"
else
    echo "File $file already exists."
fi

# create samplesheet
python3 "$file" "$fastq_dir" "$CAP_DATA_PATH/raw/samplesheet/manual_samplesheet.csv" \
    -r1 "_1.fastq.gz" \
    -r2 "_2.fastq.gz"
