#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=express
#SBATCH --time=02:00:00

# Base URLs
GTEX_ADULT="https://storage.googleapis.com/adult-gtex"
DGTEX="https://storage.googleapis.com/developmental-gtex"
BRAINSPAN_URL="https://www.brainspan.org/api/v2/well_known_file_download/267666525"

### DOWNLOAD ###
# relative paths for GTEx files
FILES=(
    "bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz"
    "annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
    "annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"
)

# relative paths for dGTEx files
DFILES=(
    "bulk-gex/v1/rna-seq/dGTEx_Analysis_2025-07-21_v1_RNASeQCv2.4.3_gene_tpm.gct.gz"
    "annotations/v1/metadata-files/dGTEx_Analysis_2026-01-30_v1_Sample_Attributes_Dataset_DS_open.xlsx"
    "annotations/v1/metadata-files/dGTEx_Analysis_2026-01-30_v1_Subject_Phenotypes_Dataset_DS_open.xlsx"
)

echo "Downloading Adult GTEx files..."
for file in "${FILES[@]}"; do
    cap_data_download --subdirectory="gtex/adult" "$GTEX_ADULT/$file"
done

echo "Downloading Developmental GTEx files..."
for dfile in "${DFILES[@]}"; do
    cap_data_download --subdirectory="gtex/developmental" "$DGTEX/$dfile"
done

# relative path for BrainSpan file
echo "Downloading BrainSpan Spatiotemporal Transcriptome..."
cap_data_download \
  --subdirectory="brainspan" \
  --source-file-name="brainspan_data.zip" \
  "$BRAINSPAN_URL"

echo "Data download complete."
echo "Unzipping BrainSpan files..."
unzip -d "${CAP_DATA_PATH}/brainspan" "${CAP_DATA_PATH}/brainspan/brainspan_data.zip"

echo "Unzipping GTEx GCT files..."
gunzip "${CAP_DATA_PATH}/gtex/"*/*.gct.gz

echo "=========================================="
echo "Download and unzipping complete!"
echo "Total runtime: $SECONDS seconds"
echo "=========================================="

