#!/bin/bash

# create sym link to reference genome directory using lab path (stored in .caprc file)
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/GENCODE_human/hg38_release_47/gencode.v47.primary_assembly.annotation.gtf"
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/GENCODE_human/hg38_release_47/GRCh38.primary_assembly.genome.fa"
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/GENCODE_mm39/release_M31/GRCm39.primary_assembly.genome.fa"
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/GENCODE_mm39/release_M31/GTF/gencode.vM31.primary_assembly.annotation.gtf"
cap_data_link "$LASSEIGNE_LAB_PATH/DATASET_dir/setbp1-hd/raw/"
