# Increase the connection buffer to 16MB
Sys.setenv("VROOM_CONNECTION_SIZE" = 32000000)
##### Permuation Analysis of SETBP1 Targets #####
## GOAL: To determine the presence and/or directional bias of SETBP1 targets in postnatal age windows.

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(ggplot2)
  library(tidyverse)
  library(vroom)
  library(data.table)
  library(patchwork)
})
###### SET UP ######
# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--indir",
  type = "character",
  help = "Path to the input data directory to access expression data"
)
parser$add_argument("-o", "--outdir",
  type = "character",
  help = "Path to the directory to write results"
)
# get command line options, if help option encountered print help and exit
args <- parser$parse_args()
indir <- args$indir
outdir <- args$outdir

# create output dir if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load helper to skip GCT headers
load_gct <- function(path, valid_sampids = NULL) {
  header <- names(fread(path, skip = 2, nrows = 0, sep = "\t", quote = "", nThread = 2))
  
  if (!is.null(valid_sampids)) {
    keep_cols <- intersect(header, c("Name", "Description", valid_sampids))
  } else {
    keep_cols <- header
  }
  
  cat("Extracting", length(keep_cols), "columns from", basename(path), "...\n")
  df <- fread(path, skip = 2, select = keep_cols, sep = "\t", quote = "")
  
  return(as_tibble(df))
}

# Function to standardize tissue names
standardize_tissue <- function(x) {
  x <- gsub(" ", "_", x)          # Spaces to underscores
  x <- gsub("-", "_", x)          # Dashes to underscores
  x <- gsub("[()]", "", x)        # Remove parentheses
  x <- gsub("__+", "_", x)        # Replace double underscores with single
  x <- trimws(x)                  # Trim whitespace
  x <- gsub("_$", "", x)          # Remove trailing underscores
  return(x)
}

# Load Metadata
cat("\nLoading adult gtex metadata...\n")
adult_meta <- vroom(paste0(indir, "gtex/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"))
adult_pheno <- vroom(paste0(indir, "gtex/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt"))

cat("\nLoading developmental gtex metadata...\n")
dev_meta   <- vroom(paste0(indir, "gtex/dGTEx_Analysis_2026-01-30_v1_Sample_Attributes_Dataset_DS_open.txt"))
dev_pheno <- vroom(paste0(indir, "gtex/dGTEx_Analysis_2026-01-30_v1_Subject_Phenotypes_Dataset_DS_open.txt"))

# SETBP1 targets
setbp1_targets_table <- read_csv(paste0(indir, "setbp1_target_pathways.csv"))
setbp1_targets <- unique(setbp1_targets_table$GeneSymbol)

# Define Phenotypically Relevant Tissues (Whitlock et al. 2024)
relevant_tissues <- c("Brain", "Muscle", "Heart", "Kidney", "Bladder", 
                      "Lung", "Small Intestine", "Stomach", "Esophagus", 
                      "Bone Marrow", "Blood")

# Extract SUBJID from SAMPID (e.g., GTEX-1117F-0226 -> GTEX-1117F)
adult_meta_clean <- adult_meta %>%
  mutate(SUBJID = sub("^([^-]+-[^-]+)-.*", "\\1", SAMPID)) %>%
  left_join(adult_pheno, by = "SUBJID") %>% # Join age/sex from phenotype file
  mutate(Tissue_Key = standardize_tissue(SMTSD)) %>%
  select(SAMPID, SUBJID, AGE, SEX, SMTS, SMTSD, Tissue_Key) %>%
  mutate(Dataset = "Adult_GTEx")

# dGTEx usually already has SUBJECT_ID, we just rename for consistency
dev_meta_clean <- dev_meta %>%
  left_join(dev_pheno, by = "SUBJECT_ID") %>%
  rename(SUBJID = SUBJECT_ID, SAMPID = SAMPLE_ID) %>%
  mutate(
    # Combine dGTEx -- some tissues have < 3 samples currently...
    AGE = "Pediatric",# case_when(
#      AGECOHORT == 1 ~ "0-2",
#      AGECOHORT == 2 ~ "2-8",
#      AGECOHORT == 3 ~ "8-13",
#      AGECOHORT == 4 ~ "13-18",
#      TRUE ~ "Unknown"
#    ),
    Tissue_Key = standardize_tissue(SMTSD),
    Dataset = "dGTEx"
  ) %>%
  select(SAMPID, SUBJID, AGE, SEX, SMTS, SMTSD, Tissue_Key, Dataset)

# Combine metadata and filter for relevant tissues
donor_census <- bind_rows(adult_meta_clean, dev_meta_clean) %>%
  filter(grepl(paste(relevant_tissues, collapse="|"), SMTS, ignore.case = TRUE)) %>%
  mutate(
SEX_Label = case_when(
SEX %in% c(1, "1") ~ "Male", 
SEX %in% c(2, "2") ~ "Female", 
TRUE ~ "Unknown")
) %>%
filter(!is.na(AGE))
cat("\nRelevant tissues found in metadata: ", unique(donor_census$SMTS))
cat("\nAge bins across data include: ", unique(donor_census$AGE))

# Specific samples to parse from TPM data
valid_samples <- donor_census$SAMPID

## Metadata Exploration ##
# Plotting the age distribution
cat("\nPlotting age distribution of gtex donors...")
p1 <- donor_census %>%
distinct(SAMPID, AGE, Dataset) %>%
  ggplot(aes(x = AGE, fill = Dataset)) +
  geom_bar(color = "white") +
  scale_y_continuous(trans = "log2") +
  theme_minimal() +
  labs(title = "Donor Distribution Across the Lifespan",
       x = "Age Group / Range",
       y = "Number of Donors") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")

print(p1)
ggsave(paste0(outdir, "gtex_metadata_donors.png"))


# Load TPM Matrices
adult_tpm <- load_gct(paste0(indir, "gtex/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct"), valid_sampids = valid_samples)
cat("\nAdult gtex count data: \n", str(adult_tpm))

dev_tpm <- load_gct(paste0(indir, "gtex/dGTEx_Analysis_2025-07-21_v1_RNASeQCv2.4.3_gene_tpm.gct"), valid_sampids = valid_samples)
cat("\nDev gtex count data: \n", str(dev_tpm))

# Clean Gene IDs (remove .version) and Join
# Keep only genes present in both datasets
cat("Filter for genes in both data sets...")

adult_tpm$Name <- sub("\\..*", "", adult_tpm$Name)
dev_tpm$Name   <- sub("\\..*", "", dev_tpm$Name)

common_genes <- intersect(adult_tpm$Name, dev_tpm$Name)
cat("\nTotal genes in common: \n", length(unique(common_genes)))

cat("\nCombining adult and dev gene counts...")
consolidated_mat <- inner_join(
  adult_tpm %>% filter(Name %in% common_genes) %>% select(-Description),
  dev_tpm %>% filter(Name %in% common_genes),
  by = "Name"
) %>% 
column_to_rownames("Name")
cat("\nCombined adult/dev gene counts structure: \n")
(str(consolidated_mat))

# Cleanup R environment
rm(adult_tpm)
rm(dev_tpm)
gc()

# Update combined_mat to use Symbols (Description)
combined_mat_symbols <- consolidated_mat %>%
  rownames_to_column("Ensembl_ID") %>%
  filter(!is.na(Description)) %>%
  # Handle potential duplicate symbols by taking the first
  distinct(Description, .keep_all = TRUE) %>% 
  column_to_rownames("Description") %>%
  select(-Ensembl_ID) %>%
  as.matrix()

# Log transform, filter zero-var genes, and row-scale
log_matrix <- log2(combined_mat_symbols + 1)
gene_sds <- apply(log_matrix, 1, sd, na.rm = TRUE)
valid_genes <- gene_sds > 0
log_matrix_filtered <- log_matrix[valid_genes, ]

# Row-scale the matrix (Mean = 0, Variance = 1 per gene)
expr_mat <- t(scale(t(log_matrix_filtered)))

# Find SETBP1 targets in combined, filtered matrix
setbp1_targets <- intersect(setbp1_targets, rownames(expr_mat))
cat("\n[Success] Found", length(setbp1_targets), "targets in the matrix.")

# Set tissue and spatio-temporal groups
analysis_groups <- donor_census %>%
  filter(SAMPID %in% colnames(expr_mat)) %>%
  group_by(SMTS, AGE) %>%
  summarise(Samples = list(SAMPID), N = n(), .groups = "drop") %>%
  filter(N >= 2) # Ensure we have at least 2 samples per test
cat("\nUnique age groups: \n", unique(donor_census$AGE))

# Permutation analysis
all_genes <- rownames(expr_mat)
n_perm <- 1000 
stats_results <- data.frame()
set.seed(123)

for(i in 1:nrow(analysis_groups)) {
  grp_tissue <- analysis_groups$SMTS[i]
  grp_age <- analysis_groups$AGE[i]
  selected_samples <- unlist(analysis_groups$Samples[i])
  
  cat("\nRunning permutations for:", grp_tissue, "| Age:", grp_age, "( N =", length(selected_samples), ")\n")
  
  select_data <- expr_mat[, selected_samples, drop = FALSE]
  
  # Observed statistic
  obs_val <- mean(rowMeans(select_data[setbp1_targets, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
  cat("\nStructure of obs values: ", str(obs_val))
  
  # Null distribution
  null_dist <- replicate(n_perm, {
    rand_genes <- sample(all_genes, length(setbp1_targets))
    mean(rowMeans(select_data[rand_genes, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
  })
  
  # Two-Sided P-values & Z-scores
  p_val_enrichment <- sum(null_dist >= obs_val) / n_perm
  p_val_depletion  <- sum(null_dist <= obs_val) / n_perm
  p_val_two_sided <- 2 * min(p_val_enrichment, p_val_depletion)
  z_score <- (obs_val - mean(null_dist)) / sd(null_dist)
  
  stats_results <- rbind(stats_results, data.frame(
    Tissue = grp_tissue,
    Age = grp_age,
    N_Samples = length(selected_samples),
    Observed = obs_val,
    Z_Score = z_score,
    P_Value = p_val_two_sided
  ))
}

# Apply FDR Correction
stats_results$FDR <- p.adjust(stats_results$P_Value, method = "fdr")
print(head(stats_results[order(stats_results$Z_Score, decreasing = TRUE), ]))

write_csv(stats_results, paste0(outdir, "postnatal_gtex_permutation_results.csv"))

## HEATMAP ##
# Order the Age column chronologically so the x-axis flows naturally
stats_results$Age <- factor(stats_results$Age, 
                            levels = c("Pediatric", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "Unknown"))

# Filter out 'Unknown' ages
plot_data <- subset(stats_results, Age != "Unknown")

# Generate the Heatmap
ggplot(plot_data, aes(x = Age, y = Tissue, fill = Z_Score)) +
  geom_tile(color = "white", size = 0.5) +
  # Add significance stars for FDR < 0.05
  geom_text(aes(label = ifelse(FDR < 0.05, "*", "")), 
            color = "black", size = 6, vjust = 0.7) +
  # Diverging color palette: Blue (depleted) -> White (0) -> Red (enriched)
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, 
                       name = "Z-Score\n(Enrichment)") +
  theme_minimal(base_size = 14) +
  labs(title = "SETBP1 Target Expression Across Tissues and Development",
       x = "Age Bracket",
       y = "Tissue") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank() # Removes grid lines for a cleaner tile look
  )

ggsave(paste0(outdir, "SETBP1_target_exp_gtex_heatmap.png"))

cat("Heatmap saved. Finished!")

