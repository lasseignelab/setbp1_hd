#!/usr/bin/env bash

# SLURM settings
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=short

# CAPTURE will download or symlink to the $CAP_DATA_PATH -- within the project's /data/ directory

#### HUMAN ####
# GRCh38 GENCODE v47 GTF annotation
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
# GRCm39 GENCODE v31 GTF annotation
cap_data_download \
  --unzip \
  --md5sum="a5efd718f75c1128d51e3e13fd85f60d" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.primary_assembly.annotation.gtf.gz" \

# Download the FASTA file
cap_data_download \
  --unzip \
  --md5sum="7079525a8ea7495e46966a64a231ae6c" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz" \
