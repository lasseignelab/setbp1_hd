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

# Loading in downloaded daa
# BrainSpan
expr_raw <- read_csv(paste0(indir, "brainspan/expression_matrix.csv"), col_names = FALSE)
rows_meta <- read_csv(paste0(indir, "brainspan/rows_metadata.csv"))
cols_meta <- read_csv(paste0(indir, "brainspan/columns_metadata.csv"))
print(paste0("Unique samples in BrainSpan: ", length(unique(cols_meta$donor_id))))
n_ages <- distinct(cols_meta, donor_id, age)
n_ages %>%
# add_count(age, donor_id) %>%
ggplot(aes(y = age)) +
geom_bar()
ggsave(paste0(outdir, "brainspan_age_donor_counts.png"))

# SETBP1 targets
setbp1_targets_table <- read_csv(paste0(indir, "setbp1_target_pathways.csv"))
setbp1_targets <- unique(setbp1_targets_table$GeneSymbol)

# Keeping track of downloaded genes
cat("Loaded", nrow(rows_meta), "genes and", nrow(cols_meta), "samples.\n")
cat("Number of SETBP1 target names read:", length(setbp1_targets), "\n\n")

# Creating Expression Matrix
row_nums_expr <- expr_raw[[1]]
if (!all(row_nums_expr == rows_meta$row_num)) {
  stop("Row numbers in expression_matrix.csv and rows_metadata.csv do not align.")
}
expr_raw_mat <- as.matrix(expr_raw[,-1])
cat("Expression matrix dimensions (genes x samples):", 
    paste(dim(expr_raw_mat), collapse = " x "), "\n")

# Transform, rm zero variance genes, and scale
log_matrix <- log2(expr_raw_mat + 1)
gene_sds <- apply(log_matrix, 1, sd)
valid_genes <- gene_sds > 0

cat("Removing", sum(!valid_genes), "genes with zero variance.\n")

log_matrix_filtered <- log_matrix[valid_genes, ]
rows_meta_filtered <- rows_meta[valid_genes, ]

expr_mat <- t(scale(t(log_matrix_filtered)))
rownames(expr_mat) <- rows_meta_filtered$gene_symbol
colnames(expr_mat) <- cols_meta$column_num
cat("Final matrix dimensions:", paste(dim(expr_mat), collapse = " x "), "\n")

# Filter SETBP1 targets for genes in brainspan data
setbp1_targets <- intersect(setbp1_targets, rownames(expr_mat))

# # Brain regions of interest for SETBP1-HD
# target_regions <- c("STR", "CBC", "VFC", "MFC", "STC")

# Converts BrainSpan age into years relative to birth (e.g., 8 pcw = 8 post-conception weeks = negative years; 4 mos = 4 months post-natal = 4/12 years; 1 yrs = 1 year after birth)
# Defining age windows and defining disease window (i.e., timeframe relevant for testing purposes); SETBP1 clinical developmental window = 0 to 2 years
cols_meta <- cols_meta %>%
  mutate(
    # Exctract numbers, then convert to numeric
    weeks  = as.numeric(str_extract(age, "\\d+")) * ifelse(grepl("pcw", age), 1, NA),
    months = as.numeric(str_extract(age, "\\d+")) * ifelse(grepl("mos", age), 1, NA),
    years  = as.numeric(str_extract(age, "\\d+")) * ifelse(grepl("yrs", age), 1, NA) 
    ) %>%
  mutate(dev_window = case_when(
           !is.na(weeks)  & weeks <= 12 ~ "Early_Prenatal",
    !is.na(weeks)  & weeks >= 13 & weeks <= 24 ~ "Mid_Prenatal",
    !is.na(weeks)  & weeks >= 25 & weeks <= 40 ~ "Late_Prenatal",
    # Increased window to capture rapid synaptogenesis and increase sample size
    (!is.na(months)) | (!is.na(years) & years <= 1) ~ "Infancy",
    !is.na(years)  & years >= 2   & years <= 5  ~ "Early_Childhood",
    !is.na(years)  & years >= 6   & years <= 16  ~ "Childhood",
    !is.na(years)  & years >= 18                ~ "Adult",
    TRUE ~ "Other" # Catches anything that doesn't fit
         )) #%>%
         # filter(structure_acronym %in% target_regions, dev_window != "Other")

cat("Sample number by age window:", table(cols_meta$dev_window))
write_csv(cols_meta, paste0(outdir, "metadata_development_windows.csv"))

cols_meta %>%
  distinct(age, donor_id, dev_window) %>%
  ggplot(aes(y = reorder(age, dev_window))) +
  geom_bar() +
  facet_wrap(~ dev_window)

ggsave(paste0(outdir, "brainspan_metadata_donors_by_age.png"))


# Permutation analyses 
## NOTE: n_perm can be changed
set.seed(123)  
n_perm <- 5000  
all_genes <- rownames(expr_mat)
gene_results <- data.frame(row.names = setbp1_targets)
stats_results <- data.frame()

# Calculating average expression for every gene for different windows
windows <- unique(cols_meta$dev_window)

# for (reg in target_regions) {
for (win in windows) {
# Samples in each window + region combo
  cat("\nRunning permutation for: \n", win)
  selected_samples <- cols_meta %>% 
      filter(
      # structure_acronym == reg, 
      dev_window == win
      ) %>% 
      pull(column_num) %>% as.character()
  
  # Ensure those samples exist in the expression matrix
    selected_samples <- intersect(selected_samples, colnames(expr_mat))

  if(length(selected_samples) < 2) next # Skip windows with insufficient data
  
  cat("\nFinding observation stat...\n")
  # Observed statistic (mean of gene set in this window)
  select_data <- expr_mat[, selected_samples, drop = FALSE]
  # Target gene means, NAs removed
  gene_means <- rowMeans(select_data[setbp1_targets, , drop=FALSE], na.rm = TRUE)
  obs_val <- mean(gene_means)
  
  cat("\nFinding null distribution...\n")
  # Null distribution
  null_dist <- replicate(n_perm, {
    rand_genes <- sample(all_genes, length(setbp1_targets))
    # Random gene means, NAs removed
    mean(rowMeans(select_data[rand_genes, , drop=FALSE], na.rm = TRUE))
  })
  
  cat("\nCalculating p-value and z-scores... \n")
  # Calculate One-Sided P-values
  p_val_enrichment <- sum(null_dist >= obs_val) / n_perm
  p_val_depletion  <- sum(null_dist <= obs_val) / n_perm
  
  # Calculate Two-Sided P-value (Recommended)
  p_val_two_sided <- 2 * min(p_val_enrichment, p_val_depletion)
  
  # Calculate Z-score
  z_score <- (obs_val - mean(null_dist)) / sd(null_dist)

  # Fix dataframe binding issue
  gene_results[[win]] <- gene_means
  
  stats_results <- rbind(stats_results, data.frame(
    window = win,
    observed = obs_val,
    z_score = z_score,
    p_val_enrichment = p_val_enrichment,
    p_val_depletion = p_val_depletion,
    p_val_two_sided = p_val_two_sided
  ))
}
# }
# Apply FDR correction for multiple hypotheses
stats_results$p_val_adj <- p.adjust(stats_results$p_val_two_sided, method = "fdr")
print(stats_results)
write_csv(stats_results, paste0(outdir, "temporal_window_perm_stats.csv"))

rownames(gene_results) <- setbp1_targets
write_csv(gene_results, paste0(outdir, "temporal_window_gene_means.csv"))
cat("\n\nPermutations complete, results saved. \n")
