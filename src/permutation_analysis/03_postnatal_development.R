# Increase the connection buffer to 16MB
Sys.setenv("VROOM_CONNECTION_SIZE" = 16000000)
##### Permuation Analysis of SETBP1 Targets #####
## GOAL: To determine the presence and/or directional bias of SETBP1 targets toward pre-natal or post-natal window.
## RESULTS: SETBP1 targets are robustly up-regulated during pre-natal window (i.e., early brain development)
## NOTES: Both mean and median were tested for these permutation analyses to determine 

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
load_gct <- function(path) {
  # skip = 2 to bypass the GCT header
  # header = TRUE to get the column names correctly
  df <- fread(path, skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
  
  # Convert to tibble if you prefer the tidyverse look
  return(as_tibble(df))
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

# --- 1. Process Adult GTEx Metadata ---
# Extract SUBJID from SAMPID (e.g., GTEX-1117F-0226 -> GTEX-1117F)
adult_meta_clean <- adult_meta %>%
  mutate(SUBJID = sub("^([^-]+-[^-]+)-.*", "\\1", SAMPID)) %>%
  left_join(adult_pheno, by = "SUBJID") %>% # Join age/sex from phenotype file
  select(SAMPID, SUBJID, AGE, SEX, SMTS, SMTSD) %>%
  mutate(Dataset = "Adult_GTEx")

# --- 2. Process Developmental GTEx Metadata ---
# dGTEx usually already has SUBJECT_ID, we just rename for consistency
# Update dev_meta_clean to decode AGECOHORT
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
    Dataset = "dGTEx"
  ) %>%
  select(SAMPID, SUBJID, AGE, SEX, SMTS, SMTSD, Dataset)

# --- 3. Combine into Master Metadata ---
master_meta <- bind_rows(adult_meta_clean, dev_meta_clean)

# --- 4. Final Cleanup: Uniform Age and Sex Labels ---
# GTEx Adult uses "1" for Male, "2" for Female; dGTEx might use strings.
donor_census <- master_meta %>%
  mutate(
    SEX_Label = case_when(
      SEX %in% c(1, "1", "Male", "M") ~ "Male",
      SEX %in% c(2, "2", "Female", "F") ~ "Female",
      TRUE ~ "Unknown"
    )
  )

# Verify the merge
table(donor_census$Dataset, donor_census$SEX_Label)

## Metadata Exploration ##
# --- Plot 1: Number of Donors per Age Group ---

# Plotting the age distribution
cat("Plotting age distribution of gtex donors...")
p1 <- ggplot(donor_census, aes(x = AGE, fill = Dataset)) +
  geom_bar(color = "white") +
  theme_minimal() +
  labs(title = "Donor Distribution Across the Lifespan",
       x = "Age Group / Range",
       y = "Number of Donors") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")

print(p1)
ggsave(paste0(outdir, "gtex_metadata_donors.png"))


# 1. Load Matrices
cat("Loading adult gtex counts...")
adult_tpm <- load_gct(paste0(indir, "gtex/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"))
cat("Loading developmental gtex counts...")
dev_tpm   <- load_gct(paste0(indir, "gtex/dGTEx_Analysis_2025-07-21_v1_RNASeQCv2.4.3_gene_tpm.gct.gz"))
cat("/nAdult gtex count data: /n", str(adult_tpm))
cat("/nDev gtex count data: /n", str(dev_tpm))

# 2. Clean Gene IDs (remove .version) and Join
# We keep only genes present in both datasets
cat("Filter for genes in both data sets...")

adult_tpm$Name <- sub("\\..*", "", adult_tpm$Name)
dev_tpm$Name   <- sub("\\..*", "", dev_tpm$Name)

common_genes <- intersect(adult_tpm$Name, dev_tpm$Name)
cat("/nTotal genes in common: /n", length(unique(common_genes)))

cat("Combining adult and dev gene counts...")
combined_mat <- inner_join(
  adult_tpm %>% filter(Name %in% common_genes) %>% select(-Description),
  dev_tpm %>% filter(Name %in% common_genes) %>% select(-Description),
  by = "Name"
) %>% 
column_to_rownames("Name") %>%
mutate(across(everything(),  ~ log2(.x + 1)))
cat("/nCombined adult/dev gene counts: /n", str(combined_mat))

# 3. Log Transformation
# combined_mat <- log2(combined_mat + 1)

# --- Plot 2: SETBP1 Target Trajectory ---
# 1. Calculate the mean expression of your targets for every sample
# Ensure combined_mat is your Log2(TPM+1) matrix
# target_means <- colMeans(combined_mat[intersect(setbp1_targets, rownames(combined_mat)), ], na.rm = TRUE)
target_means <- combined_mat %>%
filter(row.names(.) %in% setbp1_targets) %>%
rowwise() %>%
summarise(mean(c_across(everything())))
cat("\ntarget_means: \n", str(target_means))

# 2. Create a data frame for plotting
plot_data <- data.frame(
  SAMPID = names(target_means),
  Mean_Expression = target_means
) %>%
left_join(donor_census, by = "SAMPID")
#  mutate(Stage = case_when(
#    SAMPID %in% dev_meta_clean ~ "Developmental/Pediatric",
#    SAMPID %in% adult_meta_clean ~ "Adult",
#    TRUE ~ "Other"
#  )) %>%
#  filter(Stage != "Other")
cat("\nplot_data: \n", str(plot_data))

# 3. Boxplot comparison
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

plot_data_region <- plot_data %>%
   filter(SAMPID %in% region_data$SAMPID)
#  inner_join(region_data, by = "SAMPID")

ggplot(plot_data_region, aes(x = reorder(SMTSD, Mean_Expression), y = Mean_Expression, fill = AGE)) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal() +
  labs(title = "SETBP1 Target Expression by Brain Region",
       x = "",
       y = "Mean Log2(TPM+1)")

ggsave(paste0(outdir, "gtex_metadata_brain_regions.png"))

# Identify Brain Samples
# Note: GTEx uses SMTS (Tissue) and SMTSD (Detailed Tissue)
adult_brain_ids <- adult_meta_clean %>% filter(SMTS == "Brain") %>% pull(SAMPID)
dev_brain_ids   <- dev_meta_clean %>% filter(grepl("Brain", SMTS)) %>% pull(SAMPID)

# Ensure the IDs exist in our matrix
adult_cols <- intersect(adult_brain_ids, colnames(combined_mat))
dev_cols   <- intersect(dev_brain_ids, colnames(combined_mat))

cat("Running permutations...")
run_stage_permutation <- function(matrix_data, target_genes, stage_cols, n_reps = 1000) {
  
  # Subset matrix to the specific stage
  stage_mat <- as.matrix(matrix_data[, stage_cols])
  all_genes <- rownames(stage_mat)
  
  # 1. Observed Stat (Mean expression of your genes in this stage)
  # We take mean across genes, then mean across samples
  obs_stat <- mean(colMeans(stage_mat[intersect(target_genes, all_genes), , drop=FALSE]))
  
  # 2. Permutation Loop
  null_dist <- numeric(n_reps)
  set.seed(123)
  n_target <- length(intersect(target_genes, all_genes))
  
  for(i in 1:n_reps) {
    rand_genes <- sample(all_genes, n_target)
    null_dist[i] <- mean(colMeans(stage_mat[rand_genes, , drop=FALSE]))
  }
  
  p_val <- sum(null_dist >= obs_stat) / n_reps
  return(list(obs = obs_stat, p = p_val, null = null_dist))
}

# Execute for both stages
adult_results <- run_stage_permutation(combined_mat, setbp1_targets, adult_cols)
write_rds(adult_results, paste0(outdir, "gtex_adult_permutations.rds"))

dev_results   <- run_stage_permutation(combined_mat, setbp1_targets, dev_cols)
write_rds(dev_results, paste0(outdir, "devgtex_adult_permutations.rds"))

cat("Permutations complete! Results saved.")
