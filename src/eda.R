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
  source("./src/functions.R")
})

# create parser object
parser <- ArgumentParser()

# specify desired options
parser$add_argument("-r", "--rse",
  type = "character",
  metavar = "PATH",
  help = "RSE file (rds format) for EDA from recount3 data (00_data_download or 01_data_preprocessing)" # nolint
)
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

print(paste0("Reading in RSE from: ", args$rse))

#### Set variables for dataset ####
rse <- readRDS(args$rse)
study <- unique(rse@colData[["sra.study_title"]])

# save results to output dir
filepath <- paste0(args$output_dir, sra, "/")

# create dir if needed
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

#### Metadata exploration ####
meta <- as.data.frame(rse@colData) %>%
  rownames_to_column(var = "rownames") %>%
  select_if(function(col) length(unique(col)) > 1) %>%
  mutate_if(is.character, as.factor)

# frequency of metadata variables
sra_meta <- meta %>%
  select(starts_with("sra"))

pdf(paste0(filepath, "metadata_feat_freq.pdf"))
for (i in 2:ncol(sra_meta)) {
  freq <- table(sra_meta[, i])
  p <- ggplot(as.data.frame(freq), aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Freq), vjust = -.5) +
    xlab(paste(colnames(sra_meta)[i])) +
    coord_flip()
  print(p)
}
dev.off()

print(paste0(
  "Completed metadata variable frequency plots and saved: ",
  filepath, "metadata_feat_freq.pdf"
))

#### Detect/Filter Gene Counts ####

# check for counts type input and adjust counts accordingly
if (args$counts_type == "adj_cts") { # default
  if (is.null(assays(rse)$adj_cts)) {
    # auto switch to TPM if counts_type not supplied and N/A adj_cts assay
    print(paste0(
      "'counts_type' argument left as default ", args$counts_type,
      " but no 'adj_cts' assay found in RSE -- using TPM assay instead!"
    ))
    cts <- assays(rse)$TPM
    # filtering for 0 variance genes
    cts_filtered <- calc_nonzero(cts)
  } else {
    cts_filtered <- as.data.frame(assays(rse)$adj_cts)
  }
} else {
  cts <- assays(rse)[[args$counts_type]]
  # filtering for 0 variance genes
  cts_filtered <- calc_nonzero(cts)
}

#### Sample-Sample Correlation ####?
top_cts <- cts_filtered[
  c(
    order(
      apply(cts_filtered, 1, var),
      decreasing = TRUE
    )[1:5000]
  ),
]
sample_cor <- cor(top_cts)

meta_anno <- meta %>%
  column_to_rownames("rownames") %>%
  select(starts_with("sra"))

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
  )
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
  var_explained > args$var_cutoff / 100
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
  labs(title = str_wrap(paste0("Scree plot of ", sra, ": ", study)), width = 23)
# save plot
ggsave(paste0(filepath, "PC_scree_plot.png"))


# correlation of PCs to numeric QC metadata
pca_cor <- meta %>% # pc_df %>%
  # select(1:nrow(var_explained_cutoff)) %>%
  column_to_rownames(var = "rownames") %>%
  select(starts_with("recount_qc") & ends_with("both")) %>%
  # filter for data that's individual
  rename_with(~ gsub("recount_qc.star.", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("_both", "", .x, fixed = TRUE)) %>%
  bind_cols(pc_keep, .) %>%
  select(where(is.numeric)) %>%
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

pca_meta <- meta %>%
  column_to_rownames(var = "rownames") %>%
  select(contains("sra")) %>% # use sra metadata
  # filter for data that's individual
  select_if(function(col) length(unique(col)) < length(col)) %>%
  bind_cols(pc_keep, .)

# PCA color by column of interest if supplied by user
# if (!is.null(treatment)) {
if (!is.null(args$treatment)) {
  rse_pca(pca_meta, var_explained_cutoff,
    color_by = args$treatment,
    sra = "sra", study = "study",
    filepath = paste0(filepath, "PCA_coloredby_", args$treatment, ".pdf")
  )
}

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
        sra = sra,
        filepath = filepath
      )
      pca_topvar_plot(imp_feats, pca_meta, var_explained_cutoff, sra, study,
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

