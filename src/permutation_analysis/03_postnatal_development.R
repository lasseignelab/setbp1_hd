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
    # Decoding the dGTEx Age Cohorts
    AGE = case_when(
      AGECOHORT == 1 ~ "0-2",
      AGECOHORT == 2 ~ "2-8",
      AGECOHORT == 3 ~ "8-13",
      AGECOHORT == 4 ~ "13-18",
      TRUE ~ "Unknown"
    ),
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
)
cat("\nRelevant tissues found in metadata: ", unique(donor_census$SMTS))

# Verify the merge
cat("\nVerify donor census metadata merge: \n", 
table(donor_census$Dataset, donor_census$SEX_Label))

# Specific samples to parse from TPM data
valid_samples <- donor_census$SAMPID
cat("\nSamples to parse from TPM data: \n", str(valid_samples))

## Metadata Exploration ##
# Plotting the age distribution
cat("\nPlotting age distribution of gtex donors...")
p1 <- ggplot(donor_census, aes(x = AGE, fill = Dataset)) +
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

dev_tpm   <- load_gct(paste0(indir, "gtex/dGTEx_Analysis_2025-07-21_v1_RNASeQCv2.4.3_gene_tpm.gct"), valid_sampids = valid_samples)
cat("\nDev gtex count data: \n", str(dev_tpm))

# Update tissue names
colnames(adult_tpm) <- standardize_tissue(colnames(adult_tpm))

# Clean Gene IDs (remove .version) and Join
# Keep only genes present in both datasets
cat("Filter for genes in both data sets...")

adult_tpm$Name <- sub("\\..*", "", adult_tpm$Name)
dev_tpm$Name   <- sub("\\..*", "", dev_tpm$Name)

common_genes <- intersect(adult_tpm$Name, dev_tpm$Name)
cat("\nTotal genes in common: \n", length(unique(common_genes)))

cat("\nCheck adult GTEx colnames match the dGTEx tissue keys...")
print(unique(dev_meta_clean$Tissue_Key))

cat("\nCombining adult and dev gene counts...")
consolidated_mat <- inner_join(
  adult_tpm %>% filter(Name %in% common_genes) %>% select(-Description),
  dev_tpm %>% filter(Name %in% common_genes), # %>% select(-Description),
  by = "Name"
) %>% 
column_to_rownames("Name") %>%
mutate(across(where(is.numeric),  ~ log2(.x + 1)))
# mutate(across(everything(),  ~ log2(.x + 1)))
cat("\nCombined adult/dev gene counts: \n", str(consolidated_mat))

# Create a mapping table from the original data
# gene_map <- dev_tpm %>% 
#  select(Name, Description) %>%
#  distinct()

# Update combined_mat to use Symbols (Description)
combined_mat_symbols <- consolidated_mat %>%
  rownames_to_column("Ensembl_ID") %>%
#  left_join(gene_map, by = c("Ensembl_ID" = "Name")) %>%
  filter(!is.na(Description)) %>%
  # Handle potential duplicate symbols by taking the first
  distinct(Description, .keep_all = TRUE) %>% 
  column_to_rownames("Description") %>%
  select(-Ensembl_ID)

cat("\n[Fix] Matrix row names are now symbols. First few:\n")
head(rownames(combined_mat_symbols))

# Consolidate duplicate columns by averaging them
# This handles the .x and .y issue
# consolidated_mat <- combined_mat_symbols %>%
#  as.matrix()

# Rename columns to remove .x and .y
# colnames(consolidated_mat) <- sub("\\.[xy]$", "", colnames(consolidated_mat))
cat("\nCheck colnames of consolidated mat: \n", head(colnames(consolidated_mat)))

# Average columns with the same name
consolidated_mat <- sapply(unique(colnames(consolidated_mat)), function(col) {
  if (sum(colnames(consolidated_mat) == col) > 1) {
    rowMeans(consolidated_mat[, colnames(consolidated_mat) == col, drop=FALSE])
  } else {
    consolidated_mat[, colnames(consolidated_mat) == col]
  }
})

cat("\n[Fix] Consolidated columns. Total unique tissues:", ncol(consolidated_mat))

# Calculate the mean expression of your targets for every sample
target_genes_present <- intersect(setbp1_targets, rownames(combined_mat_symbols))
cat("\n[Success] Found", length(target_genes_present), "targets in the matrix.")

# Calculate means: result is a named numeric vector
target_means <- colMeans(consolidated_mat[target_genes_present, , drop = FALSE], na.rm = TRUE)
cat("\ntarget_means: \n", str(target_means))

# Create a data frame for plotting
plot_data <- data.frame(
  Tissue_Key = names(target_means),
  Mean_Expression = as.numeric(target_means)
) %>%
left_join(
donor_census %>% 
select(Tissue_Key, Dataset, AGE, SMTS, SMTSD) %>% 
distinct(), by = "Tissue_Key"
)
cat("\nplot_data: \n", str(plot_data))

# Boxplot comparison
p2 <- ggplot(plot_data, aes(x = AGE, y = Mean_Expression, fill = Dataset)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_classic() +
  labs(title = "SETBP1 Target Expression: Dev vs. Adult",
       subtitle = "Integrated GTEx + dGTEx Brain Samples",
       y = "Mean Log2(TPM+1)") +
  scale_fill_manual(values = c("Adult" = "#2c3e50", "Developmental/Pediatric" = "#e74c3c"))

print(p2)
ggsave(paste0(outdir, "gtex_stage_exp_boxplot.png"))

# Map samples to specific brain regions (SMTSD in GTEx)
region_data <- bind_rows(
  adult_meta_clean %>% select(SAMPID, SMTSD),
  dev_meta_clean %>% select(SAMPID, SMTSD)
)
cat("Plotting expression by brain region...")

# Filter for brain tissues
plot_data_brain <- plot_data %>%
  filter(SMTS == "Brain") %>%
  filter(!is.na(Mean_Expression) & !is.nan(Mean_Expression))

ggplot(
plot_data_brain, aes(
x = reorder(SMTSD, Mean_Expression), y = Mean_Expression, fill = AGE)
) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal() +
  labs(title = "SETBP1 Target Expression by Brain Region",
       subtitle = "Median Expression per Tissue/Age Group",
       x = "Brain Region",
       y = "Mean Log2(TPM+1)")

ggsave(paste0(outdir, "gtex_metadata_brain_regions.png"))

# Identify Brain Samples
# SMTS = Tissue and SMTSD = Detailed Tissue
adult_brain_ids <- adult_meta_clean %>% 
filter(SMTS == "Brain") %>% 
pull(SAMPID)
dev_brain_ids   <- dev_meta_clean %>% 
filter(grepl("Brain", SMTS)) %>% 
pull(SAMPID)

# Ensure the IDs exist in our matrix
adult_cols <- intersect(adult_brain_ids, colnames(combined_mat))
dev_cols   <- intersect(dev_brain_ids, colnames(combined_mat))

cat("Running permutations...")
run_median_permutation <- function(matrix_data, target_genes, tissue_colnames, n_reps = 1000) {
  
# Filter matrix for the specific tissues (e.g., all columns starting with "Brain")
  # Use drop=FALSE to maintain matrix structure even if only 1 tissue matches
  sub_mat <- as.matrix(matrix_data[, tissue_colnames, drop = FALSE])
  all_genes <- rownames(sub_mat)
  
  # Observed Statistic: 
  valid_targets <- intersect(target_genes, all_genes)
  n_target <- length(valid_targets)
  
  if(n_target == 0) return(list(obs = NA, p = NA, message = "No genes found"))
  
  obs_stat <- mean(rowMeans(sub_mat[valid_targets, , drop = FALSE]))
  
  # Permutation Loop
  null_dist <- numeric(n_reps)
  set.seed(123)
  
  for(i in 1:n_reps) {
    rand_genes <- sample(all_genes, n_target)
    null_dist[i] <- mean(rowMeans(sub_mat[rand_genes, , drop = FALSE]))
  }
  
  # Calculate P-value
  p_val <- sum(null_dist >= obs_stat) / n_reps
  
  return(list(obs = obs_stat, p = p_val, null = null_dist, n = n_target))
}

# Identify Brain columns for Adult and Dev
brain_cols <- colnames(consolidated_mat)[grepl("Brain", colnames(consolidated_mat))]

# Ensure these columns actually exist in the combined matrix
adult_brain_cols <- intersect(brain_cols, adult_meta_clean$Tissue_Key)
dev_brain_cols   <- intersect(brain_cols, dev_meta_clean$Tissue_Key)

cat("\n[Permutation] Testing", length(adult_brain_cols), "Adult brain tissues and", 
    length(dev_brain_cols), "Dev brain tissues.")

# Run Tests
adult_res <- run_median_permutation(consolidated_mat, setbp1_targets, adult_brain_cols)
dev_res <- run_median_permutation(consolidated_mat, setbp1_targets, dev_brain_cols)

# Clean Output Summary
cat("\n==========================================")
cat("\nPERMUTATION TEST RESULTS (SETBP1 Targets)")
cat("\n==========================================")
cat("\nADULT BRAIN:")
cat("\n  Observed Mean Exp:", round(adult_res$obs, 4))
cat("\n  Empirical P-value:", adult_res$p)
cat("\n  Significant?     :", ifelse(adult_res$p < 0.05, "YES", "NO"))

cat("\n\nDEVELOPMENTAL BRAIN:")
cat("\n  Observed Mean Exp:", round(dev_res$obs, 4))
cat("\n  Empirical P-value:", dev_res$p)
cat("\n  Significant?     :", ifelse(dev_res$p < 0.05, "YES", "NO"))
cat("\n==========================================\n")


# Execute for both stages
# adult_results <- run_stage_permutation(combined_mat, setbp1_targets, adult_cols)
write_rds(adult_res, paste0(outdir, "gtex_adult_permutations.rds"))

# dev_results   <- run_stage_permutation(combined_mat, setbp1_targets, dev_cols)
write_rds(dev_res, paste0(outdir, "devgtex_adult_permutations.rds"))

# cat("Permutations complete! Results saved.")
