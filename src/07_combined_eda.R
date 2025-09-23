suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(gprofiler2)
  source(paste0(getwd(), "/src/functions.R"))
})
# create parser object
parser <- ArgumentParser()

# specify desired options
parser$add_argument("-f", "--functions",
  type = "character",
  metavar = "PATH",
  default = NULL,
  help = "Path to functions file to source custom functions" # nolint
)
parser$add_argument("-r", "--rse",
  type = "character",
  metavar = "PATH",
  help = "Filepath to directory with *_prepped_rse.rds files"
)
parser$add_argument("-v", "--var_cutoff",
  type = "integer",
  default = 3,
  help = "Integer value representing % variance threshold cutoff for PC's to retain (i.e. 3 = retain PC's explaining > 3.0% total variation)" # nolint
)
parser$add_argument("-o", "--output_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory to create and write EDA results"
)

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# source functions file if supplied
if (!is.null(args$functions)) source(args$functions)

# start timer
ptm <- proc.time()

#### Set variables for dataset ####
print(paste0("Reading in RSE from: ", args$rse))

var_cutoff <- args$var_cutoff
sra <- args$sra
filepath <- args$output_dir
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

c_rse <- readRDS(paste0(args$rse, "cardo_prepped_rse.rds"))
s_rse <- readRDS(paste0(args$rse, "shaw_prepped_rse.rds"))
t1_rse <- readRDS(paste0(args$rse, "tanaka1_prepped_rse.rds"))
t2_rse <- readRDS(paste0(args$rse, "tanaka2_prepped_rse.rds"))

meta_comb <- t1_rse@metadata %>%
  bind_rows(c_rse@metadata, s_rse@metadata) %>%
  mutate(author = ifelse(
    BioProject == "PRJNA746989", "Cardo",
    ifelse(
      BioProject == "PRJNA1093043", "Shaw", "Tanaka"
    )
  ))

write_csv(meta_comb, file = paste0(filepath, "combined_metadata.csv"))

# combine cardo and shaw data, subset for nonzero var genes
nonzero_genes <- intersect(rownames(c_rse), rownames(s_rse))
cts_comb <- cbind(c_rse[nonzero_genes, ], s_rse[nonzero_genes, ])
cts_filtered <- assay(cts_comb, "salmon.merged.gene_tpm")

# convert tanaka mouse IDs to human orthologs
conv <- gorth(
  rownames(t1_rse),
  source_organism = "mmusculus",
  target_organism = "hsapiens", filter_na = TRUE, mthreshold = 1
)

# subset RSEs by genes w/ mouse-human ortholog available
ortho_genes <- intersect(rownames(cts_filtered), conv$ortholog_ensg)
ortho_filtered <- cts_filtered[rownames(cts_filtered) %in% ortho_genes, ] %>%
  rownames_to_column()

# mouse genes w/ multiple orthologs for 1 human
comb_genes <- t1_rse@assays@data@listData$salmon.merged.gene_tpm %>%
  rownames_to_column() %>%
  mutate(rowname = gsub("\\.\\d+$", "", rowname)) %>%
  inner_join(conv, by = c("rowname" = "input")) %>%
  inner_join(ortho_filtered, by = c("ortholog_ensg" = "rowname"))

comb_cts <- comb_genes %>%
  select(ortholog_ensg | starts_with("SRR"))

# top var human genes
top_cts <- comb_cts[
  c(
    order(
      apply(comb_cts, 1, var),
      decreasing = TRUE
    )[1:1000]
  ),
]


top_cts_filt <- top_cts %>%
  # remove non-unique genes
  distinct(ortholog_ensg, .keep_all = TRUE) %>%
  column_to_rownames("ortholog_ensg")

### SAVE
write_csv(top_cts_filt, file = paste0(filepath, "combined_counts_topvar.csv"))

# calculate PC's
cts_mat <- t(top_cts_filt)
cts_mat <- cts_mat[, apply(cts_mat, 2, var) > 0]
# pc <- prcomp(t(top_cts_filt), scale = TRUE)
pc <- prcomp(cts_mat, scale = TRUE)
pc_df <- as.data.frame(pc$x)

# calculate variance explained
var_explained_df <- data.frame(
  PC = paste0("PC", seq_along(pc_df)),
  var_explained = pc$sdev^2 / sum(pc$sdev^2)
)

# select variance cutoff to filter by for retaining PC's
var_explained_cutoff <- filter(
  var_explained_df,
  var_explained > 2 / 100
)
# filter PC's for cutoff
pc_keep <- pc_df %>%
  select(seq_len(nrow(var_explained_cutoff)))
print(paste0("Summary of top var-gene PCs across data sets: "))
print(summary(pc_keep))

pca_meta <- meta_comb %>%
  filter(rownames(.) %in% rownames(pc_keep)) %>%
  filter(rownames(.) %in% colnames(top_cts_filt)) %>%
  # filter metadata for features of interest
  select(BioProject, cells, genotype, Instrument, rep, condition, author) %>%
  bind_cols(pc_keep, .) %>%
  mutate(
    Organism = ifelse(
      BioProject %in% c("PRJNA898600", "PRJNA974068"), "Mouse", "Human"
    ),
    across(where(is.character), as.factor)
  )

meta_feats <- setdiff(colnames(pca_meta), colnames(pc_keep))
print(paste0("Metadata features to color by:"))
meta_feats

# color by condition
for (i in seq_along(meta_feats)) {
  rse_pca(pca_meta,
    var_explained_cutoff,
    color_by = meta_feats[i],
    shape = meta_feats[6],
    sra = "Combined Datasets",
    filepath = paste0(
      filepath, "PCA_coloredby_", meta_feats[i], "_WTvsHD.pdf"
    )
  )
}

# color by cells
for (i in seq_along(meta_feats)) {
  rse_pca(pca_meta,
    var_explained_cutoff,
    color_by = meta_feats[i],
    shape = meta_feats[2],
    sra = "Combined Datasets",
    filepath = paste0(
      filepath, "PCA_coloredby_", meta_feats[i], "_cells.pdf"
    )
  )
}

# color by author
for (i in seq_along(meta_feats)) {
  rse_pca(pca_meta,
    var_explained_cutoff,
    color_by = meta_feats[i],
    shape = meta_feats[7],
    sra = "Combined Datasets",
    filepath = paste0(
      filepath, "PCA_coloredby_", meta_feats[i], "_author.pdf"
    )
  )
}
