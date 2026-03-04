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
expr_mat <- as.matrix(expr_raw[,-1])
cat("Expression matrix dimensions (genes x samples):", 
    paste(dim(expr_mat), collapse = " x "), "\n")
rownames(expr_mat) <- rows_meta$gene_symbol
colnames(expr_mat) <- cols_meta$column_num
cat(str(expr_mat))

write_csv(as.data.frame(expr_mat), paste0(outdir, "aggregate_gene_exp_mat.csv"))

# Converts BrainSpan age into years relative to birth (e.g., 8 pcw = 8 post-conception weeks = negative years; 4 mos = 4 months post-natal = 4/12 years; 1 yrs = 1 year after birth)
parse_age_to_years <- function(age_str) {
  age_str <- trimws(age_str)
  if (grepl("pcw", age_str)) {
    w <- as.numeric(sub(" pcw", "", age_str))
    return((w - 40) / 52)
  } else if (grepl("mos", age_str)) {
    m <- as.numeric(sub(" mos", "", age_str))
    return(m / 12)
  } else if (grepl("yrs", age_str)) {
    y <- as.numeric(sub(" yrs", "", age_str))
    return(y)
  } else {
    stop(paste("Unknown age format:", age_str))
  }
}

# Defining age windows and defining disease window (i.e., timeframe relevant for testing purposes); SETBP1 clinical developmental window = 0 to 2 years
cols_meta <- cols_meta %>%
  mutate(
    age_years = vapply(age, parse_age_to_years, numeric(1)),
    prenatal  = age_years < 0,
    postnatal = age_years >= 0,
    disease_window = age_years >= 0 & age_years <= 2
  )

cat("Prenatal samples:", sum(cols_meta$prenatal),
    " | Postnatal samples:", sum(cols_meta$postnatal),
    " | Disease-window samples:", sum(cols_meta$disease_window), "\n\n")


# Calculating average expression for every gene for different windows
## Indexing samples by age window
prenatal_idx  <- cols_meta$prenatal
postnatal_idx <- cols_meta$postnatal
disease_idx   <- cols_meta$disease_window
if (!any(prenatal_idx) | !any(postnatal_idx)) {
  stop("No prenatal or no postnatal samples detected – check age parsing.")
}
## Row means for each gene in each window
prenatal_means <- rowMeans(expr_mat[, prenatal_idx, drop = FALSE], na.rm = TRUE)
postnatal_means <- rowMeans(expr_mat[, postnatal_idx, drop = FALSE], na.rm = TRUE)

## Temporal stat to determine whether gene expression is higher or lower prenatally
prenatal_minus_postnatal <- prenatal_means - postnatal_means

## Creating table to hold all info
gene_stats <- tibble(
  gene_symbol = rownames(expr_mat),
  prenatal_mean = prenatal_means,
  postnatal_mean = postnatal_means,
  prenatal_minus_postnatal = prenatal_minus_postnatal#,
  # disease_mean = disease_means
)
## Filter out genes with all-NAs or degenerate values
gene_stats <- gene_stats %>%
  filter(!is.na(prenatal_minus_postnatal))

cat("Genes with valid prenatal/postnatal means:", nrow(gene_stats), "\n\n")
write_csv(gene_stats, paste0(outdir, "prenatal_postnatal_gene_means.csv"))

# Creating master list of all of the genes and filtering for SETBP1 targets
master_stats <- gene_stats$gene_symbol
setbp1_stats <- intersect(setbp1_targets, master_stats)
cat("SETBP1 targets found in BrainSpan universe:", 
    length(setbp1_stats), "of", length(setbp1_targets), "\n\n")

if (length(setbp1_stats) < 5) {
  warning("Very few SETBP1 targets found in BrainSpan; results may be unstable.")
}
## Stats vector for future use
stat_vec <- gene_stats$prenatal_minus_postnatal
names(stat_vec) <- gene_stats$gene_symbol

## Creating observed Stat for SETBP1 targets for permutation analyses
obs_stat <- mean(stat_vec[setbp1_stats], na.rm = TRUE)
cat("Observed mean (prenatal - postnatal) for SETBP1 targets:", obs_stat, "\n\n")

# Permutation analyses 
## NOTE: n_perm can be changed
set.seed(123)  
n_perm <- 5000  
set_size <- length(setbp1_stats)
all_genes <- setdiff(gene_stats$gene_symbol, setbp1_targets)
perm_stats <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  rand_genes <- sample(all_genes, size = set_size, replace = FALSE)
  perm_stats[i] <- mean(stat_vec[rand_genes], na.rm = TRUE)
  if (i %% 1000 == 0) {
    cat("Permutation", i, "of", n_perm, "complete\r")
  }
}
cat("\n\nPermutation loop finished.\n")

## Empirical p-value to test whether SETBP1 targets more upregulated during the prenatal window instead of post natal
p_emp <- (sum(perm_stats >= obs_stat) + 1) / (n_perm + 1)
null_mean <- mean(perm_stats)
null_sd   <- sd(perm_stats)
z_score   <- (obs_stat - null_mean) / null_sd

cat("=== Permutation results (prenatal - postnatal) ===\n")
cat("Observed mean (SETBP1):", obs_stat, "\n")
cat("Null mean:", null_mean, "| Null SD:", null_sd, "\n")
cat("Z-score:", z_score, "\n")
cat("Empirical p (SETBP1 > random):", p_emp, "\n\n")

# Plot for testing purposes
png(paste0("histogram_mean_perm.png"), width = 20, height = 15, units = "cm", res = 300)
hist(perm_stats, breaks = 50,
      main = "Null distribution of mean (prenatal - postnatal)",
     xlab = "Mean (prenatal - postnatal) over random gene sets")
abline(v = obs_stat, col = "red", lwd = 2)
legend("topright", legend = c("Observed SETBP1 mean"), lwd = 2, col = "red")
dev.off()


# Median based permutation tests to confirm SETBP1 targets, as a whole, display pre-natal upregulated gene expression

# Data Cleaning
## Medians of each row, separated by prenatal/postnatal windo
prenatal_medians  <- apply(expr_mat[, prenatal_idx, drop = FALSE], 1, median, na.rm = TRUE)
postnatal_medians <- apply(expr_mat[, postnatal_idx, drop = FALSE], 1, median, na.rm = TRUE)
prenatal_minus_postnatal_med <- prenatal_medians - postnatal_medians

## Gene stats median table 
gene_stats_med <- tibble(
  gene_symbol = rownames(expr_mat),
  diff_median = prenatal_minus_postnatal_med
) %>% filter(!is.na(diff_median))
write_csv(gene_stats_med, paste0(outdir, "prenatal_postnatal_gene_median.csv"))

## Median Vector List
median_stats_vector <- gene_stats_med$diff_median
names(median_stats_vector) <- gene_stats_med$gene_symbol

## Median Observed Stat
median_obs_stat <- median(median_stats_vector[setbp1_stats], na.rm = TRUE)

# Permutation test for median
set.seed(123)
perm_stats_med <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  # Pick random genes
  rand_genes <- sample(names(median_stats_vector), size = set_size, replace = FALSE)
  perm_stats_med[i] <- median(median_stats_vector[rand_genes], na.rm = TRUE)
}

p_emp_med <- (sum(perm_stats_med >= median_obs_stat) + 1) / (n_perm + 1)
z_score_med <- (median_obs_stat - mean(perm_stats_med)) / sd(perm_stats_med)

# Print the results
cat("=== Median-based Results ===\n")
cat("Observed Median:", median_obs_stat, "\n")
cat("Empirical p-value:", p_emp_med, "\n")
cat("Z-score:", z_score_med, "\n")


# Check how many unique values are in your random trials
cat("Unique values in perm_stats_med:", length(unique(perm_stats_med)), "\n")

# Plot the Median results
## NOTE: Red line may not show up on graph due to x-axis stopping at 0
png(paste0("histogram_median_perm.png"), width = 20, height = 15, units = "cm", res = 300)
hist(perm_stats_med, breaks = 50,
     main = "Median-based Null Distribution",
     xlab = "Median (prenatal - postnatal) of random sets",
     col = "lightblue")
abline(v = median_obs_stat, col = "red", lwd = 2)
legend("topright", legend = c("Observed SETBP1 Median"), lwd = 2, col = "red")
dev.off()


