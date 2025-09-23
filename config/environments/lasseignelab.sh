#!/bin/bash

# create sym link to reference genome directory using lab path (stored in .caprc file)
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/GENCODE_human/hg38_release_45/gencode.v45.annotation.gtf"
cap_data_link "$LASSEIGNE_LAB_PATH/GENOME_dir/GENCODE_human/hg38_release_45/GRCh38.primary_assembly.genome.fa"

# nf-core RNA-seq outputs
cap_data_link "$LASSEIGNE_LAB_PATH/PROJECT_dir/setbp1_hd/cardo"
cap_data_link "$LASSEIGNE_LAB_PATH/PROJECT_dir/setbp1_hd/shaw"
cap_data_link "$LASSEIGNE_LAB_PATH/PROJECT_dir/setbp1_hd/tanaka"
