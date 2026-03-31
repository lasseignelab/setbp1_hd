# Increase the connection buffer to 16MB
Sys.setenv("VROOM_CONNECTION_SIZE" = 32000000)
##### Permuation Analysis of SETBP1 Targets #####
## GOAL: To determine the presence and/or directional bias of
## SETBP1 targets across tissues and ages.

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(ggplot2)
  library(tidyverse)
  library(vroom)
  library(data.table)
  library(readxl)
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
addir <- paste0(indir, "gtex/adult/GTEx_Analysis_")
devdir <- paste0(indir, "gtex/developmental/dGTEx_Analysis_")

outdir <- args$outdir

# create output dir if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load helper to skip GCT headers
load_gct <- function(path, valid_sampids = NULL) {
  header <- names(
    fread(path, skip = 2, nrows = 0, sep = "\t", quote = "", nThread = 2)
  )

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
  x <- gsub(" ", "_", x) # Spaces to underscores
  x <- gsub("-", "_", x) # Dashes to underscores
  x <- gsub("[()]", "", x) # Remove parentheses
  x <- gsub("__+", "_", x) # Replace double underscores with single
  x <- trimws(x) # Trim whitespace
  x <- gsub("_$", "", x) # Remove trailing underscores
}

# Load Metadata
cat("\nLoading adult gtex metadata...\n")
adult_meta <- vroom(
  paste0(addir, "v11_Annotations_SampleAttributesDS.txt")
)
adult_pheno <- vroom(
  paste0(addir, "v11_Annotations_SubjectPhenotypesDS.txt")
)

cat("\nLoading developmental gtex metadata...\n")
dev_meta <- read_xlsx(
  paste0(
    devdir,
    "2026-01-30_v1_Sample_Attributes_Dataset_DS_open.xlsx"
  )
)
dev_pheno <- read_xlsx(
  paste0(
    devdir,
    "2026-01-30_v1_Subject_Phenotypes_Dataset_DS_open.xlsx"
  )
)

# SETBP1 targets
setbp1_targets_table <- read_csv(paste0(indir, "260303_setbp1_targets.csv"))
setbp1_targets <- unique(setbp1_targets_table$Target)

# Define Phenotypically Relevant Tissues (Whitlock et al. 2024)
relevant_tissues <- c(
  "Brain", "Muscle", "Heart", "Kidney", "Bladder",
  "Lung", "Small Intestine", "Stomach", "Esophagus",
  "Bone Marrow", "Blood"
)

# Extract SUBJID from SAMPID (e.g., GTEX-1117F-0226 -> GTEX-1117F)
adult_meta_clean <- adult_meta %>%
  mutate(SUBJID = sub("^([^-]+-[^-]+)-.*", "\\1", SAMPID)) %>%
  left_join(adult_pheno, by = "SUBJID") %>% # Join age/sex from phenotype file
  mutate(Tissue_Key = standardize_tissue(SMTSD)) %>%
  select(SAMPID, SUBJID, AGE, SEX, SMTS, SMTSD, Tissue_Key) %>%
  mutate(Dataset = "Adult_GTEx")

# Rename for consistency
dev_meta_clean <- dev_meta %>%
  left_join(dev_pheno, by = "SUBJECT_ID") %>%
  rename(SUBJID = SUBJECT_ID, SAMPID = SAMPLE_ID) %>%
  mutate(
    # Combine dGTEx -- some tissues have < 3 samples currently...
    AGE = "Pediatric",
    Tissue_Key = standardize_tissue(SMTSD),
    Dataset = "dGTEx"
  ) %>%
  select(SAMPID, SUBJID, AGE, SEX, SMTS, SMTSD, Tissue_Key, Dataset)

# Combine metadata and filter for relevant tissues
donor_census <- bind_rows(adult_meta_clean, dev_meta_clean) %>%
  filter(
    grepl(
      paste(relevant_tissues, collapse = "|"), SMTS,
      ignore.case = TRUE
    )
  ) %>%
  mutate(
    SEX_Label = case_when(
      SEX %in% c(1, "1") ~ "Male",
      SEX %in% c(2, "2") ~ "Female",
      TRUE ~ "Unknown"
    )
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
  labs(
    title = "Donor Distribution Across the Lifespan",
    x = "Age Group / Range",
    y = "Number of Donors"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")

print(p1)
ggsave(paste0(outdir, "gtex_metadata_donors.png"))


# Load TPM Matrices
adult_tpm <- load_gct(
  paste0(
    addir, "2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct"
  ),
  valid_sampids = valid_samples
)

dev_tpm <- load_gct(
  paste0(
    devdir, "2025-07-21_v1_RNASeQCv2.4.3_gene_tpm.gct"
  ),
  valid_sampids = valid_samples
)

# Clean Gene IDs (remove .version), join, and keep only genes in both datasets
cat("Filter for genes in both data sets...")

adult_tpm$Name <- sub("\\..*", "", adult_tpm$Name)
dev_tpm$Name <- sub("\\..*", "", dev_tpm$Name)

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

### SETBP1 TOP CO-EXPRESSED TARGETS ###
# Subset for SETBP1 and target gene expression
setbp1_expr <- as.numeric(expr_mat["SETBP1", ])
target_mat <- expr_mat[setbp1_targets, ]

# Calculate Spearman correlation between SETBP1 and all targets
results_mat <- apply(target_mat, 1, function(target_expr) {
  test <- cor.test(setbp1_expr, target_expr, method = "spearman", exact = FALSE)
  c(Spearman_Rho = test$estimate[["rho"]], P_Value = test$p.value)
})

# Data wrangling and FDR correction
coexpression_df <- as.data.frame(t(results_mat)) %>%
  rownames_to_column("Target_Gene") %>%
  mutate(FDR = p.adjust(P_Value, method = "fdr")) %>%
  arrange(desc(Spearman_Rho))

# View top correlated targets
cat("Most Positively Correlated:\n")
print(head(coexpression_df, 10))
cat("\nMost Negatively Correlated:\n")
print(tail(coexpression_df, 10))

write_csv(coexpression_df, paste0(outdir, "gtex_target_coexpression.csv"))

### SEX-BIASED EXPRESSION ###
# Calculate the mean expression of SETBP1 targets for every single sample
sample_target_means <- colMeans(
  expr_mat[setbp1_targets, , drop = FALSE],
  na.rm = TRUE
)

# Filter for known sexes
sex_data <- donor_census %>%
  filter(
    SAMPID %in% names(sample_target_means),
    SEX_Label %in% c("Male", "Female")
  ) %>%
  mutate(Target_Mean = sample_target_means[SAMPID])

write_csv(sex_data, paste0(outdir, "gtex_target_sex_expression.csv"))
cat("\nTarget expression by sex saved to: gtex_target_sex_expression.csv")

# Run the Wilcoxon rank sum test grouped by Tissue and Age
sex_results <- sex_data %>%
  group_by(SMTS, AGE) %>%
  # Ensure the group has at least 1 male and 1 female to compare
  filter(sum(SEX_Label == "Male") >= 3 & sum(SEX_Label == "Female") >= 3) %>%
  summarise(
    Male_N = sum(SEX_Label == "Male"),
    Female_N = sum(SEX_Label == "Female"),
    P_Value = wilcox.test(Target_Mean ~ SEX_Label, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  # Apply FDR correction
  mutate(FDR = p.adjust(P_Value, method = "fdr")) %>%
  arrange(P_Value)

write_csv(sex_results, paste0(outdir, "gtex_sex_wilcox_results.csv"))
cat("\nSex comparison results saved to: postnatal_gtex_sex_results.csv")


### PERMUTATION ANALYSIS ACROSS TISSUES ###
# Set tissue and spatio-temporal groups
analysis_groups <- donor_census %>%
  filter(SAMPID %in% colnames(expr_mat)) %>%
  group_by(SMTS, AGE) %>%
  summarise(Samples = list(SAMPID), N = n(), .groups = "drop") %>%
  filter(N >= 2) # Ensure we have at least 2 samples per test
cat("\nUnique age groups: \n", unique(donor_census$AGE))

# Permutation analysis
all_genes <- rownames(expr_mat)
n_perm <- 10000
stats_results <- data.frame()
set.seed(123)

for (i in seq_len(nrow(analysis_groups))) {
  grp_tissue <- analysis_groups$SMTS[i]
  grp_age <- analysis_groups$AGE[i]
  selected_samples <- unlist(analysis_groups$Samples[i])

  cat(
    "\nRunning permutations for:", grp_tissue,
    "| Age:", grp_age, "( N =", length(selected_samples), ")\n"
  )

  select_data <- expr_mat[, selected_samples, drop = FALSE]

  # Observed statistic
  obs_val <- mean(
    rowMeans(select_data[setbp1_targets, , drop = FALSE], na.rm = TRUE),
    na.rm = TRUE
  )

  # Null distribution
  null_dist <- replicate(n_perm, {
    rand_genes <- sample(all_genes, length(setbp1_targets))
    mean(
      rowMeans(
        select_data[rand_genes, , drop = FALSE],
        na.rm = TRUE
      ),
      na.rm = TRUE
    )
  })

  # Two-Sided P-values & Z-scores
  p_val_enrichment <- sum(null_dist >= obs_val) / n_perm
  p_val_depletion <- sum(null_dist <= obs_val) / n_perm
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

write_csv(
  stats_results, paste0(outdir, "gtex_permutation_results.csv")
)
cat("Finished!\nResults saved to: gtex_permutation_results.csv")

#### END ####

# end timer
fptm <- proc.time()
cat("\n\nTotal runtime:  ", (fptm[3] / 60), "\n")

#  session info
cat("\n\nSession Info:\n")
print(sessionInfo())
