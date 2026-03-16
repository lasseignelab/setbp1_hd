#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=express
#SBATCH --time=02:00:00

# Base URLs
GTEX_ADULT="https://storage.googleapis.com/adult-gtex"
DGTEX="https://storage.googleapis.com/developmental-gtex"
BRAINSPAN_URL="https://www.brainspan.org/api/v2/well_known_file_download/267666525"

### Download Adult GTEx and dGTEx ###
# relative paths for GTEx files
FILES=(
    "bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
    "bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz"
    "annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
    "annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"
)

# relative paths for dGTEx files
DFILES=(
    "bulk-gex/v1/rna-seq/dGTEx_Analysis_2025-07-21_v1_RNASeQCv2.4.3_gene_tpm.gct.gz"
    "annotations/v1/metadata-files/dGTEx_Analysis_2026-01-30_v1_Sample_Attributes_Dataset_DS_open.txt"
    "annotations/v1/metadata-files/dGTEx_Analysis_2026-01-30_v1_Subject_Phenotypes_Dataset_DS_open.txt"
)

echo "Downloading Adult GTEx files..."
for file in "${FILES[@]}"; do
    cap_data_download --subdirectory="gtex" "$GTEX_ADULT/$file"
done

echo "Downloading Developmental GTEx files..."
for dfile in "${DFILES[@]}"; do
    cap_data_download --subdirectory="gtex" "$DGTEX/$dfile"
done

### Download BrainSpan Data ###
echo "Downloading BrainSpan Spatiotemporal Transcriptome..."
cap_data_download \
  --subdirectory="brainspan" \
  --source-file-name="brainspan_data.zip" \
  "$BRAINSPAN_URL"

echo "Data download complete. Unzipping brainspan files..."
mkdir -p "${CAP_DATA_PATH}/brainspan" 
unzip -d "${CAP_DATA_PATH}/brainspan" "${CAP_DATA_PATH}/brainspan/brainspan_data.zip"
echo "Unzipping GTEx files..."
gunzip "${CAP_DATA_PATH}/gtex/"*.gct.gz
# Adult GTEx
# median gene-level TPM by tissue
# bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz
# bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz
# annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt
# annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt

# Developmental GTEx TPMs and metadata
# annotations/v1/metadata-files/dGTEx_Analysis_2026-01-30_v1_Sample_Attributes_Dataset_DS_open.txt
# annotations/v1/metadata-files/dGTEx_Analysis_2026-01-30_v1_Subject_Phenotypes_Dataset_DS_open.txt

# Brainspan data
# https://www.brainspan.org/api/v2/well_known_file_download/267666525

