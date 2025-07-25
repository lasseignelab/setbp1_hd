#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=medium
#SBATCH --array=0-2

module load parallel
cap_container -c singularity "ncbi/sra-tools:3.2.1"

# print
echo "$CAP_DATA_PATH"


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
echo "$out_dir"

#### READ ACCESSION LIST ####
 while IFS= read -r accession; do
    echo "Processing $accession for $dataset ..."

    # check if files already exist in directory
    if [ -z "$(find "$out_dir" -type f -name "$accession_file*.fastq.gz")" ]; then
    # fetch fastqs for each SRR in accession list
        singularity exec \
          --bind "$out_dir":/output_dir \
          --bind "$CAP_DATA_PATH":/tmp_dir \
    	    "$CAP_CONTAINER_PATH"/sra-tools_3.2.1.sif \
    	    sh -c "prefetch --output-directory /tmp_dir/ --max-size 50g $accession && \
                fasterq-dump $accession --threads $SLURM_CPUS_PER_TASK -t /tmp_dir -O /output_dir/ && \
                rm -r /tmp_dir/$accession && \
                cd /output_dir/"
    else
        echo "Directory $data_dir already exists and has fastq files in it." >&2
    fi
done < "$accession_file"

#### COMPRESS FASTQ FILES AFTER ALL DOWNLOADS ####
echo "All files fetched for $dataset. Compressing fastq files now ..."
parallel -j "$SLURM_CPUS_PER_TASK" gzip -n ::: "$out_dir"/*.fastq
echo "All .fastq files have been compressed for $dataset."
