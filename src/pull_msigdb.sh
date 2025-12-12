#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=medium
#SBATCH --partition=express

cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=WP_PI3KAKT_SIGNALING&fileType=TSV"
cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=BIOCARTA_TGFB_PATHWAY&fileType=TSV"
cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=BIOCARTA_SHH_PATHWAY&fileType=TSV"
cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=BIOCARTA_IGF1R_PATHWAY&fileType=TSV"
cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=BIOCARTA_EGF_PATHWAY&fileType=TSV"
cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_MYC_TARGETS_V1&fileType=TSV"
cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_MYC_TARGETS_V2&fileType=TSV"
cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_TGF_BETA_SIGNALING&fileType=TSV"
cap_data_download \
  --subdirectory "mSigDB" \
  "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_WNT_BETA_CATENIN_SIGNALING&fileType=TSV"

