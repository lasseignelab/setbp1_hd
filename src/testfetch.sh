#!/bin/bash

#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=medium

cap_container -c singularity "ncbi/sra-tools:3.2.1"

#### DATASET CONFIGURATION ####
declare -A datasets=(
    ["cardo"]="$CAP_PROJECT_PATH/doc/Cardo_SRR_Acc_List.txt"
    ["shaw"]="$CAP_PROJECT_PATH/doc/Shaw_SRR_Acc_List.txt"
    ["tanaka"]="$CAP_PROJECT_PATH/doc/Tanaka_SRR_Acc_List.txt"
)

dataset_names=("cardo" "shaw" "tanaka")
dataset=${dataset_names[$SLURM_ARRAY_TASK_ID]}
accession_file=${datasets[$dataset]}

#### OUTPUT DIRECTORY ####
data_dir="$CAP_DATA_PATH/raw/$dataset"
mkdir -p "$data_dir"

#### READ ACCESSION LIST ####
while IFS= read -r SRR_ID; do
    echo "Downloading $SRR_ID for dataset $dataset..."

    # Check if files already exist
    if [ -z "$(find "$data_dir" -type f -name "$SRR_ID*.fastq.gz")" ]; then
        # Format sample number (e.g., S01, S02...)
        sample_number=$(printf "S%02d" $((SLURM_ARRAY_TASK_ID + 1)))

        singularity exec \
            --bind "$data_dir":/output_dir \
            --bind "$CAP_DATA_PATH":/tmp_dir \
            "$CAP_CONTAINER_PATH"/sra-tools_3.2.1.sif \
            sh -c "prefetch --output-directory /tmp_dir/ --max-size 50g $SRR_ID && \
                   fasterq-dump $SRR_ID -t /tmp_dir/ -O /output_dir/ && \
                   rm -r /tmp_dir/$SRR_ID && \
                   cd /output_dir/ && \
                   mv \"${SRR_ID}_1.fastq\" \"${SRR_ID}_${sample_number}_R1_001.fastq\" && \
                   mv \"${SRR_ID}_2.fastq\" \"${SRR_ID}_${sample_number}_R2_001.fastq\""
    else
        echo "FASTQ files for $SRR_ID already exist in $data_dir. Skipping." >&2
    fi
done < "$accession_file"

#### COMPRESS FASTQ FILES AFTER ALL DOWNLOADS ####
echo "All downloads complete. Compressing FASTQ files..."
module load parallel
parallel -j "$SLURM_CPUS_PER_TASK" gzip ::: "$data_dir"/*.fastq
echo "Compression complete."

