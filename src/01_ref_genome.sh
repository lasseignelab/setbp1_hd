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
  --md5sum="81cf2337bb31aeaddba7059a3b0f14f5" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz" \

# Download the FASTA files
cap_data_download \
  --unzip \
  --md5sum="0c33152bea6f09b89fb5e3f89558771e" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz" \

#### MOUSE ####
# GENCODE v39 GTF and FASTA from 10x
cap_data_download \
  --unzip \
  --md5sum="37c51137ccaeabd4d151f80dc86ce0b3" \
  "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
