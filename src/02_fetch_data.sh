#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=24G
#SBATCH --partition=medium


#### CONDA ENVIRONMENT ####
ml Anaconda3/2023.07-2

conda activate sra-env
if ! conda activate sra-env; then
    echo "Failed to activate Conda environment 'sra-env'"
    exit 1
fi
echo "Activated Conda environment 'sra-env'"


#### OUTPUT DIRECTORY ####
out_dir="$CAP_DATA_PATH/raw/fastq"
mkdir -p "$out_dir"

#### FETCH FASTQ FILES ####
# fetch fastqs for each SRR in accession list
while IFS= read -r accession; do
    echo "Processing $accession ..."
    fasterq-dump "$accession" --outdir "$out_dir"
    vdb-validate "$accession"
    echo "SRA $accession validated ..."

done < "$CAP_PROJECT_PATH/doc/SRR_Acc_List.txt"

echo "All files fetched. Compressing fastq files now ..."

#### GZIP FASTQ FILES ####
for file in "$out_dir"/*.fastq; do
    if [ -f "$file" ]; then
        gzip "$file"
        echo "Compressed $file"
    fi
done

echo "All .fastq files have been compressed."
