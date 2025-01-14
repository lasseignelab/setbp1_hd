#!/usr/bin/env bash

# SLURM settings
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=short

# CAPTURE will download to the $CAP_DATA_PATH -- default environment uses project's /data/ directory
# GENCODE v45 GTF annotation
cap_data_download \
    --md5sum="b6eeb6c9791b7a43a5504a654ff09d9a" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz"

# Primary assembly FASTA
cap_data_download \
    --md5sum="ad62ff4d71d0b5b8d7feabbec5ce86bf" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
