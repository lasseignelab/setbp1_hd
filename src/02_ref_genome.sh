#!/usr/bin/env bash

# SLURM settings
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=short

# CAPTURE will download or symlink to the $CAP_DATA_PATH -- within the project's /data/ directory

#### HUMAN ####
# GENCODE v45 GTF annotation
cap_data_download \
    --unzip \
    --md5sum="b6eeb6c9791b7a43a5504a654ff09d9a" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz"

# Primary assembly FASTA
cap_data_download \
    --unzip \
    --md5sum="ad62ff4d71d0b5b8d7feabbec5ce86bf" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"

#### MOUSE ####
# GENCODE v39 GTF and FASTA from 10x
cap_data_download \
  --unzip \
  --md5sum="37c51137ccaeabd4d151f80dc86ce0b3" \
  "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
