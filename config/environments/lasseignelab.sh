#!/bin/bash

# create sym link to reference genome directory using lab path (stored in .caprc file)
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/GENCODE_human/hg38_release_45/gencode.v45.annotation.gtf"
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/GENCODE_human/hg38_release_45/GRCh38.primary_assembly.genome.fa"
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/refdata-gex-GRCm39-2024-A"
