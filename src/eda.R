suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(tidymodels)
  library(DALEX)
  library(ranger)
  library(gbm)
  library(e1071)
  source(paste0(getwd(), "/functions.R"))
})

# create parser object
parser <- ArgumentParser()

# specify desired options
parser$add_argument("-r", "--rse",
  type = "character",
  metavar = "PATH",
  help = "RSE file (rds format) for EDA from recount3 data (00_data_download or 01_data_preprocessing)" # nolint
)
parser$add_argument("-s", "--sra",
                    type = "character",
                    metavar = "PATH",
                    help = "SRA or dataset for naming purposes" # nolint
)
# parser$add_argument("-m", "--meta",
#                     type = "character",
#                     metavar = "PATH",
#                     help = "Metadata file in csv format" # nolint
# )
parser$add_argument("-v", "--var_cutoff",
  type = "integer",
  default = 3,
  help = "Integer value representing % variance threshold cutoff for PC's to retain (i.e. 3 = retain PC's explaining > 3.0% total variation)" # nolint
)
parser$add_argument("-t", "--treatment",
  type = "character",
  default = NULL,
  help = "Treatment/condition of interest in RSE colData col to color PCA plots by" # nolint
)
parser$add_argument("-o", "--output_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory to create and write EDA results"
)
# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# start timer
ptm <- proc.time()

#### Set variables for dataset ####
print(paste0("Reading in RSE from: ", args$rse))

rse <- readRDS(args$rse)
meta <- rse@metadata
# meta <- read.csv(args$meta, row.names = 1, stringsAsFactors = TRUE)
var_cutoff <- args$var_cutoff
sra <- args$sra
# meta <- read.csv("./doc/tanaka_metadata.csv", row.names = 1, stringsAsFactors = TRUE)
# rse <- readRDS("./data/tanaka/star_salmon/salmon.merged.gene_counts.SummarizedExperiment.rds")
# meta <- rse@metadata
# filepath <- "./results/archive/"
# var_cutoff <- 4
# sra <- "tanaka"

# create output dir if needed
filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

#### Metadata exploration ####
# frequency of metadata variables

groupings <- meta %>%
  # filter out metadata that's completely unique
  select_if(function(col) length(unique(col)) < length(col)) %>%
  # filter out data that's completely redundant
  select_if(function(col) length(unique(col)) != 1) 

pdf(paste0(filepath, "metadata_feat_freq.pdf"))
for (i in 2:ncol(groupings)) {
  freq <- as.data.frame(table(groupings[, i]))
  p <- freq %>%
    ggplot(aes(x = freq[,1], y = Freq)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Freq), vjust = -.5) +
    xlab(paste(colnames(groupings)[i])) +
    coord_flip()
  print(p)
}
dev.off()

print(paste0(
  "Completed metadata variable frequency plots and saved: ",
  filepath, "metadata_feat_freq.pdf"
))

# #### Detect/Filter Gene Counts ####
# filtering for 0 variance genes
cts_filtered <- calc_nonzero(rse@assays@data@listData[["salmon.merged.gene_tpm"]])

#### Sample Variation ####
# cts_filtered2 <- cts_filtered[1:100,]

cts_long <- cts_filtered %>%
  rownames_to_column(var = "genes") %>%
  pivot_longer(cols = !genes, names_to = "sample_ID", values_to = "cts")

meta_cts_long <- groupings %>%
  rownames_to_column(var = "sample_ID") %>%
  left_join(cts_long, by = "sample_ID")

meta_cts_long %>%
  mutate(rep = factor(rep)) %>%
  ggplot(aes(log2(cts+1), colour = rep)) +
  geom_freqpoly(binwidth = 5) +
  facet_grid(cells ~ genotype)

# save plot
ggsave(paste0(filepath, sra, "_sample_cts_variation_.png"), height = 7, width = 14, units = "cm")

#### Sample-Sample Correlation ####
top_cts <- cts_filtered[
  c(
    order(
      apply(cts_filtered, 1, var),
      decreasing = TRUE
    )[1:5000]
  ),
]
sample_cor <- cor(top_cts)

meta_anno <- groupings[match(colnames(sample_cor), row.names(groupings)),]

png(
  paste0(
    filepath, "QC_sample_to_sample_corr.png"
  ),
  width = 35, height = 35, units = "cm", res = 300
)
ComplexHeatmap::Heatmap(as.matrix(sample_cor),
  top_annotation = HeatmapAnnotation(
    df = meta_anno,
    show_legend = FALSE
  ),
  cluster_rows = FALSE, cluster_columns = FALSE
)
dev.off()

print(paste0(
  "Completed sample-sample correlations by top 5000 var genes and plotted: ",
  filepath, "QC_sample_to_sample_corr.png"
))

#### PCA ####

# calculate PC's
pc <- prcomp(t(cts_filtered), scale = TRUE)
pc_df <- as.data.frame(pc$x)

# calculate variance explained
var_explained_df <- data.frame(
  PC = paste0("PC", seq_along(pc_df)),
  var_explained = pc$sdev^2 / sum(pc$sdev^2)
)
summary(var_explained_df)

# select variance cutoff to filter by for retaining PC's
var_explained_cutoff <- filter(
  var_explained_df,
  var_explained > var_cutoff / 100
)

# if threshold only retains 1 PC -- automatically retain PC1 AND PC2
if (nrow(var_explained_cutoff) < 2) {
  var_explained_cutoff <- filter(
    var_explained_df,
    PC %in% c("PC1", "PC2")
  )
}

# filter PC's for cutoff
pc_keep <- pc_df %>%
  select(1:nrow(var_explained_cutoff))

# Scree plot
var_explained_cutoff %>%
  ggplot(aes(
    x = reorder(PC, var_explained, decreasing = TRUE), y = var_explained
  )) +
  geom_point(size = 4) +
  # geom_line()+
  xlab("PC's Above Var Explained Threshold") +
  labs(title = str_wrap(paste0("Scree plot of ", ": ")), width = 23)
# save plot
ggsave(paste0(filepath, "PC_scree_plot.png"))


# correlation of PCs to numeric QC metadata
pca_cor <- meta %>% # pc_df %>%
  # select(1:nrow(var_explained_cutoff)) %>%
  select_if(is.numeric) %>%
  select_if(function(col) length(unique(col)) != 1) %>%
  # select(starts_with("recount_qc") & ends_with("both")) %>%
  # filter for data that's individual
  # rename_with(~ gsub("recount_qc.star.", "", .x, fixed = TRUE)) %>%
  # rename_with(~ gsub("_both", "", .x, fixed = TRUE)) %>%
  bind_cols(pc_keep, .) %>%
  # select(where(is.numeric)) %>%
  cor()


png(
  paste0(
    filepath, "QC_PC_correlation.png"
  ),
  width = 30, height = 30, units = "cm", res = 300
)
ComplexHeatmap::Heatmap(as.matrix(pca_cor))
dev.off()

print(paste0(
  "Completed PC-recount QC correlations and plotted: ",
  filepath, "QC_PC_correlation.png"
))

#### Metadata features driving variation ####

pca_meta <- groupings %>%
  # filter for data that's individual
  # select_if(function(col) length(unique(col)) < length(col)) %>%
  bind_cols(pc_keep, .) %>%
  mutate(across(where(is.character), as.factor))

# PCA color by column of interest(s) if supplied by user
if (sra == "cardo"){
  int_feats <- c("genotype", "cells", "differentiation_days", "sample", "sample_type", "condition")
}
if (sra == "shaw"){
  int_feats <- c("genotype", "cells", "sample", "sample_type", "ACMG", "condition")
}
if (sra == "tanaka"){
  int_feats <- c("genotype", "AGE", "sample", "sample_type", "sex", "condition")
}

for (i in seq_along(int_feats)){
  rse_pca(pca_meta, var_explained_cutoff,
          color_by = int_feats[i],
          sra = sra,
          filepath = paste0(filepath, "PCA_coloredby_", int_feats[i], ".pdf")
  )
}
# rse_pca(pca_meta, var_explained_cutoff,
#         color_by = "ACMG",
#         sra = "sra",
#         filepath = paste0(filepath, "PCA_coloredby_sampletype.pdf")
# )


#### Feature importance of each PC ####

# model types to use for predicting feature importance for each PC
model_types <- c("rf", "svm", "gbm")

for (i in seq_along(model_types)) {
  tryCatch(
    {
      print(paste0("Creating ", model_types[i], " and explainer"))
      imp_feats <- get_imp_feats(pca_meta,
        modeltype = model_types[i],
        suffix = "PC",
        numpc = nrow(var_explained_cutoff),
        sra = "test",
        filepath = filepath
      )
      pca_topvar_plot(imp_feats, pca_meta, var_explained_cutoff, sra,
        modeltype = model_types[i], filepath
      )
      print(paste0(
        "PCAs colored by most-important-features saved: ",
        filepath, model_types[i], "_feature_importance.pdf"
      ))
    },
    error = function(e) {
      cat("ERROR ", conditionMessage(e), "\n")
    }
  )
}


#### Compare feature importance by model ####

comp_mods <- list()
for (i in seq_along(model_types)) {
  tryCatch(
    {
      comp_mods[[i]] <- readRDS(
        paste0(
          filepath, "feature_importance_", model_types[i], ".rds"
        )
      )
    },
    error = function(e) {
      cat(
        "Model comparison could not be plotted: 1 or more failed models\n",
        conditionMessage(e), "\n"
      ) # no lint
    }
  )
}
if (length(comp_mods) == 3) {
  pdf(paste0(filepath, "PC_feat_imp_model_comparison.pdf"))
  for (i in seq_along(var_explained_cutoff)) {
    print(plot(comp_mods[[1]][[i]], comp_mods[[2]][[i]], comp_mods[[3]][[i]]))
  }
  dev.off()
}

# end timer
fptm <- proc.time()
(fptm[3] / 60)

