# SETBP1-HD: Intersection of Regulatory Analysis and Signature Reversion Uncovers Therapeutic Drugs and Targets
## Purpose
This repository contains the code and analysis pipeline for identifying repurposable drugs for * *SETBP1* * Haploinsufficiency Disorder (* *SETBP1* *-HD). By deriving a consensus transcriptomic signature from * *in-vitro* * models and integrating it with * *SETBP1* * regulatory targets, we prioritize FDA-approved therapeutics (such as celecoxib and buspirone) capable of reversing disease-associated molecular disruptions.

## Dependencies
This codebase was built using the [CAPTURE](https://github.com/lasseignelab/capture) framework (cap version 1.0) and uses a combination of Conda environments and Docker containers to ensure reproducibility. Below are the environments/images required and the scripts that run within them:

- **Conda Environment (bin/conda/nfcore-env.yml)**
  - 04_human_preprocess.sh
  - 04_mouse_preprocess.sh
- **Docker: ncbi/sra-tools:3.2.1**
  - 02_fetch_data.sh
- **Docker: lizzyr/lw_models:0.1.1**
  - 05_data_cleaning.R
  - 06_dataset_eda.R
  - 07_combined_eda.R
- **Docker: lizzyr/lw_dea:0.3.1**
  - 08_find_degs.R
  - 09_analyze_degs.R
  - 10_consensus_signature.R
  - permutation_analysis/02_job_permutation.sh
  - permutation_analysis/03_job_permutation.sh
  - permutation_analysis/04_job_prioritize.sh
- **Docker: lizzyr/sigsearch:1.0.0**
  - 11_signature_reversion.R
  - 12_drugtarget_degs.R
- **Docker: projectassistant/sitemap-scraper:latest**
  - 13_drug_safety.py

## Data
- **Signature reversion drug repurposing analyses:**
  - GSE262710: **Cardo et al.** [Identifying SETBP1 haploinsufficiency molecular pathways to improve patient diagnosis using induced pluripotent stem cells and neural disease modelling](https://pmc.ncbi.nlm.nih.gov/articles/PMC11443744/)
  - GSE180185: **Shaw et al.** [Impaired neurogenesis and neural progenitor fate choice in a human stem cell model of SETBP1 disorder](https://pmc.ncbi.nlm.nih.gov/articles/PMC9940404/)
  - PRJNA898600 (Lin−c-Kit+ cells) & PRJNA974068 (Lin−c-Kit+ Sca-1+ cells): **Tanaka et al.** [SETBP1 is dispensable for normal and malignant hematopoiesis](https://www.nature.com/articles/s41375-023-01970-5)
  - GSE301238 (fibroblasts) & GSE301239 (induced neurons): **Wong et al.** [SETBP1 variants outside the degron disrupt DNA-binding, transcription and neuronal differentiation capacity to cause a heterogeneous neurodevelopmental disorder](https://www.nature.com/articles/s41467-025-64074-x)

- **Regulatory target prioritization & permutation testing:**
   - [GTEx](https://gtexportal.org/home/downloads/adult-gtex/overview) 2025-08-22 v11 RNASeQCv2.4.3 TPM and metadata files [v11 Sample and Subject Attributes DS]
   - [dGTEx](https://gtexportal.org/home/downloads/developmental-gtex) 2025-07-21 v1 RNASeQCv2.4.3 TPM and metadata files [v1 2026-01-30 Sample and Subject Attributes])
   - [BrainSpan](https://www.brainspan.org/static/download.html) api v2, developmental transcriptome zip file with expression matrix, row metadata [genes], column metadata [samples])


## Scripts
The src/ directory contains the complete pipeline for data fetching, preprocessing, exploratory data analysis (EDA), differential expression analysis (DEA), and drug prioritization.
```
src/
├── 01_ref_genome.sh                # Downloads and builds reference genome indices.
├── 02_fetch_data.sh                # Retrieves raw RNA-seq data from public repositories.
├── 03_fastq_samplesheet_job.sh     # Generates sample sheets mapping the downloaded FASTQ files.
├── 04_human_preprocess.sh          # Preprocesses raw human reads (QC, alignment, quantification).
├── 04_mouse_preprocess.sh          # Preprocesses raw mouse reads (QC, alignment, quantification).
├── 05_data_cleaning.R              # Cleans and formats raw count matrices for downstream analysis.
├── 05_job_cleaning.sh
├── 06_dataset_eda.R                # Performs Exploratory Data Analysis (EDA) on individual datasets.
├── 06_job_eda.sh
├── 07_combined_eda.R               # Evaluates cross-study batch effects and combined data structure.
├── 07_job_combined_eda.sh
├── 08_find_degs.R                  # Conducts differential expression analysis (DEA) to find significant gene changes.
├── 08_job_dea.sh
├── 09_analyze_degs.R               # Analyzes and visualizes the DEGs (e.g., volcano plots, functional enrichment).
├── 09_job_analyze_degs.sh
├── 10_consensus_signature.R        # Derives the robust consensus transcriptomic signature across SETBP1-HD models.
├── 10_job_consensus_sig.sh
├── 11_signature_reversion.R        # Queries the LINCS database to find drugs that reverse the consensus signature.
├── 11_job_sigsearch.sh
├── 12_drugtarget_degs.R            # Intersects drug candidate targets with prioritized SETBP1 regulatory targets.
├── 12_job_drugtarget_degs.sh
├── 13_drug_safety.py               # Scrapes and evaluates pediatric/pregnancy safety profiles for top drug candidates.
├── 13_job_drug_safety.sh
├── figures/                        # Scripts for generating final publication figures.
├── functions.R                     # Custom reusable R functions used throughout the pipeline.
└── permutation_analysis/
    ├── 01_fetch_temporal_exp.sh      # Fetches BrainSpan and GTEx temporal gene expression data.
    ├── 02_job_permutation.sh
    ├── 02_spatiotemporal_perm.R      # Performs permutation testing to identify critical developmental windows for SETBP1 targets.
    ├── 03_job_permutation.sh
    ├── 03_postnatal_development.R    # Analyzes SETBP1 target expression shifts during the prenatal-to-postnatal transition.
    ├── 04_job_prioritize.sh
    └── 04_prioritize_targets.R       # Prioritizes targets based on spatiotemporal co-expression, core pathways, and known druggability.
```

## Authors 

Elizabeth J. Wilk, Sasha Taluri, Tabea M. Soelter, and Brittany N. Lasseigne


## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<!-- markdownlint-disable-next-line MD013 -->

<img src="https://static.wixstatic.com/media/29dcd0_945b269e399141df9792e4ce78f832ff~mv2.png" width="75" height="75">

## Funding

We are grateful to our funders for this work: the [SETBP1 Society](https://www.setbp1.org/) (2024 Microgrant) and the [UAB Pilot Center for Precision Animal Modeling (C-PAM)](https://sites.uab.edu/cpam/) (U54-OD030167)

## Acknowledgements

We thank the Lasseigne Lab members for their thoughtful input and discussion. We also extend our deep gratitude to the * *SETBP1* * Society and the broader community of researchers, clinicians, patients, and caregivers. Your dedication, open sharing of data, and collaborative discussions have been instrumental in advancing this research and improving the care of those affected by * *SETBP1* * variation.

## License

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.
