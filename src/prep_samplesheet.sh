#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --partition=express


cd "$CAP_DATA_PATH/raw/samplesheet/"

# select just first 4 columns for new samplesheet
cut -d ',' -f1-4 nonrna_rm_samplesheet.csv > all_samplesheet.csv

# append with manual samplesheet, removing header and adding quotes
tail -n +2 manual_samplesheet.csv | awk -F, '{print "\"" $1 "\",\"" $2 "\",\"" $3 "\",\"" $4 "\""}' >> all_samplesheet.csv
