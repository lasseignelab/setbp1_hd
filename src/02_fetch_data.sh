#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=24G
#SBATCH --partition=medium
#SBATCH --array=0-2

#### CONDA ENVIRONMENT ####
ml Anaconda3/2023.07-2

conda activate sra-env
if ! conda activate sra-env; then
    echo "Failed to activate Conda environment 'sra-env'"
    exit 1
fi
echo "Activated Conda environment 'sra-env'"

#### ACCESSION ID FILES AND OUTPUT DIRECTORIES ####
declare -A datasets=(
    ["cardo"]="$CAP_PROJECT_PATH/doc/Cardo_SRR_Acc_List.txt"
    ["shaw"]="$CAP_PROJECT_PATH/doc/Shaw_SRR_Acc_List.txt"
    ["tanaka"]="$CAP_PROJECT_PATH/doc/Tanaka_SRR_Acc_List.txt"
)

dataset_names=("cardo" "shaw" "tanaka")
dataset=${dataset_names[$SLURM_ARRAY_TASK_ID]}
accession_file=${datasets[$dataset]}

#### OUTPUT DIRECTORY ####
out_dir="$CAP_DATA_PATH/raw/$dataset"
mkdir -p "$out_dir"

#### FETCH FASTQ FILES ####
# fetch fastqs for each SRR in accession list
while IFS= read -r accession; do
    echo "Processing $accession for $dataset ..."
    fasterq-dump "$accession" --outdir "$out_dir"
    vdb-validate "$accession"
    echo "SRA $accession validated ..."

done < "$accession_file"

echo "All files fetched for $dataset. Compressing fastq files now ..."

#### GZIP FASTQ FILES ####
for file in "$out_dir"/*.fastq; do
    if [ -f "$file" ]; then
        gzip "$file"
        echo "Compressed $file for $dataset"
    fi
done

echo "All .fastq files have been compressed for $dataset."
