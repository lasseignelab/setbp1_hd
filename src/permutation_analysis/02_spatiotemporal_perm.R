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

# GTEx (rows = genes, cols = samples, col 1 = row_num)
# GCT files have 2 lines of header metadata before the actual table
load_gct <- function(path) {
  vroom(path, skip = 2, col_types = c(Name = "c", Description = "c", .default = "d"))
}

# Example: Loading Developmental GTEx TPMs
# gtex_tpm <- load_gct(paste0(indir, "gtex/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct"))
dgtex_tpm <- read_delim(paste0(indir, "gtex/dGTEx_Analysis_2026-01-30_v1_Sample_Attributes_Dataset_DS_open.txt"))
print(paste0("Unique samples in dGTEx: ", length(unique(dgtex_tpm$SUBJECT_ID))))

# 3. Log2 Transformation for your Permutation Test
# We add a pseudocount of 1 to avoid log(0)
# gtex_log_tpm <- gtex_tpm %>%
#   column_to_rownames("Name") %>%
#   select(-Description) %>%
#   mutate(across(everything(), ~ log2(.x + 1)))
# gtex_raw_expression <- read_csv("expression_matrix.csv", col_names = FALSE)

# SETBP1 targets
setbp1_targets_table <- read_csv(paste0(indir, "260303_setbp1_targets.csv"))
setbp1_targets <- unique(setbp1_targets_table$Target)
# setbp1_targets   <- read_excel("250609_setbp1_targets_names_only.xlsx", col_names = TRUE)
# setbp1_targets <- setbp1_targets[[1]] |> unique() |> as.character()

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
# Transform and scale
log_matrix <- log2(expr_raw_mat + 1)
expr_mat <- t(scale(t(log_matrix)))
rownames(expr_mat) <- rows_meta$gene_symbol
colnames(expr_mat) <- cols_meta$column_num

# Filter SETBP1 targets for genes in brainspan data
setbp1_targets <- intersect(setbp1_targets, rownames(expr_mat))

# Brain regions of interest for SETBP1-HD
target_regions <- c("STR", "CBC", "VFC", "MFC", "STC")

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
           !is.na(weeks)  & weeks >= 8  & weeks <= 13 ~ "Early_Prenatal",
    !is.na(weeks)  & weeks >= 16 & weeks <= 19 ~ "Mid_Prenatal",
    !is.na(weeks)  & weeks >= 21 & weeks <= 37 ~ "Late_Prenatal",
    # Increased window to capture rapid synaptogenesis and increase sample size
    !is.na(months) | !is.na(years) & years <= 1 ~ "Infancy",
    # Increased window to capture synaptic pruning and increase sample size
    !is.na(years)  & years >= 2   & years <= 16  ~ "Childhood",
    !is.na(years)  & years >= 18                ~ "Adult",
    TRUE ~ "Other" # Catches anything that doesn't fit
         )) %>%
         filter(structure_acronym %in% target_regions, dev_window != "Other")

cat("Sample number by age window:", table(cols_meta$dev_window))
write_csv(cols_meta, paste0(outdir, "metadata_development_windows.csv"))

# Permutation analyses 
## NOTE: n_perm can be changed
set.seed(123)  
n_perm <- 5000  
all_genes <- rownames(expr_mat)
stats_results <- data.frame()

# Calculating average expression for every gene for different windows
windows <- unique(cols_meta$dev_window)

for (reg in target_regions) {
for (win in windows) {
# Samples in each window + region combo
  cat("\nRunning permutation for: \n", reg, win)
  selected_samples <- cols_meta %>% 
      filter(structure_acronym == reg, dev_window == win) %>% 
      pull(column_num) %>% as.character()
  
  # Ensure those samples exist in the expression matrix
    selected_samples <- intersect(selected_samples, colnames(expr_mat))
  if(length(selected_samples) < 2) next # Skip windows with insufficient data
  
  cat("\nFinding observation stat...\n")
  # Observed statistic (mean of gene set in this window)
  select_data <- expr_mat[, selected_samples, drop = FALSE]
  # Target gene means, NAs removed
  obs_val <- mean(rowMeans(select_data[setbp1_targets, , drop=FALSE], na.rm = TRUE))
  
  cat("\nFinding null distribution...\n")
  # Null distribution
  null_dist <- replicate(n_perm, {
    rand_genes <- sample(all_genes, length(setbp1_targets))
    # Random gene means, NAs removed
    mean(rowMeans(select_data[rand_genes, , drop=FALSE], na.rm = TRUE))
  })
  
  cat("\nCalculating p-value and z-scores... \n")
  # Calculate P-value and Enrichment Score (Z-score)
  p_val <- sum(null_dist >= obs_val) / n_perm
  z_score <- (obs_val - mean(null_dist)) / sd(null_dist)
  
  stats_results <- rbind(stats_results, data.frame(
    region = reg,
    window = win,
    observed = obs_val,
    p_value = p_val,
    z_score = z_score
  ))
}
}
write_csv(stats_results, paste0(outdir, "temporal_window_perm_stats.csv"))
cat("\n\nPermutation loop finished.\n")

# Plot for testing purposes
# Ensure windows are in chronological order for the plot
stats_results$window <- factor(stats_results$window, levels = windows)

ggplot(stats_results, aes(x = window, y = z_score, group = region, color = region)) +
  geom_line(size = 1) +
  geom_point(aes(shape = p_value < 0.05), size = 3.5) +
  facet_wrap(~region) +
  scale_shape_manual(values = c(1, 19), name = "p < 0.05") +
  theme_minimal() +
  labs(title = "Regional Trajectories of Gene Set Enrichment",
       y = "Enrichment Z-Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(outdir, "temporal_exp_trajectory.png"))
cat("Plot saved to temporal_exp_trajectory.png")

# # Median based permutation tests to confirm SETBP1 targets, as a whole, display pre-natal upregulated gene expression
#
# # Data Cleaning
# ## Medians of each row, separated by prenatal/postnatal windo
# prenatal_medians  <- apply(expr_mat[, prenatal_idx, drop = FALSE], 1, median, na.rm = TRUE)
# postnatal_medians <- apply(expr_mat[, postnatal_idx, drop = FALSE], 1, median, na.rm = TRUE)
# prenatal_minus_postnatal_med <- prenatal_medians - postnatal_medians
#
# ## Gene stats median table 
# gene_stats_med <- tibble(
#   gene_symbol = rownames(expr_mat),
#   diff_median = prenatal_minus_postnatal_med
# ) %>% filter(!is.na(diff_median))
# write_csv(gene_stats_med, paste0(outdir, "prenatal_postnatal_gene_median.csv"))
#
# ## Median Vector List
# median_stats_vector <- gene_stats_med$diff_median
# names(median_stats_vector) <- gene_stats_med$gene_symbol
#
# ## Median Observed Stat
# median_obs_stat <- median(median_stats_vector[setbp1_stats], na.rm = TRUE)
#
# # Permutation test for median
# set.seed(123)
# perm_stats_med <- numeric(n_perm)
#
# for (i in seq_len(n_perm)) {
#   # Pick random genes
#   rand_genes <- sample(names(median_stats_vector), size = set_size, replace = FALSE)
#   perm_stats_med[i] <- median(median_stats_vector[rand_genes], na.rm = TRUE)
# }
#
# p_emp_med <- (sum(perm_stats_med >= median_obs_stat) + 1) / (n_perm + 1)
# z_score_med <- (median_obs_stat - mean(perm_stats_med)) / sd(perm_stats_med)
#
# # Print the results
# cat("=== Median-based Results ===\n")
# cat("Observed Median:", median_obs_stat, "\n")
# cat("Empirical p-value:", p_emp_med, "\n")
# cat("Z-score:", z_score_med, "\n")
#
#
# # Check how many unique values are in your random trials
# cat("Unique values in perm_stats_med:", length(unique(perm_stats_med)), "\n")
#
# # Plot the Median results
# ## NOTE: Red line may not show up on graph due to x-axis stopping at 0
# png(paste0("histogram_median_perm.png"), width = 20, height = 15, units = "cm", res = 300)
# hist(perm_stats_med, breaks = 50,
#      main = "Median-based Null Distribution",
#      xlab = "Median (prenatal - postnatal) of random sets",
#      col = "lightblue")
# abline(v = median_obs_stat, col = "red", lwd = 2)
# legend("topright", legend = c("Observed SETBP1 Median"), lwd = 2, col = "red")
# dev.off()
#
#
