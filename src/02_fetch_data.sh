#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=medium
#SBATCH --array=0-2

module load parallel
cap_container -c singularity "ncbi/sra-tools:3.2.1"

#### SRA TOOLKIT CONFIG ####
export SRA_CACHE_DIR="$CAP_DATA_PATH/tmp"
mkdir -p "$SRA_CACHE_DIR"

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
    singularity exec \
	--bind "$CAP_DATA_PATH" \
	"$CAP_CONTAINER_PATH"/sra-tools_3.2.1.sif \
    	sh -c "fasterq-dump --threads $SLURM_CPUS_PER_TASK -t $SRA_CACHE_DIR/tmp -x $accession --outdir $out_dir
    	echo \"Validating ...\"
    	vdb-validate $accession"

done < "$accession_file"

echo "All files fetched for $dataset. Compressing fastq files now ..."

#### GZIP FASTQ FILES ####
parallel -j "$SLURM_CPUS_PER_TASK" gzip ::: "$out_dir"/*.fastq

echo "All .fastq files have been compressed for $dataset."
