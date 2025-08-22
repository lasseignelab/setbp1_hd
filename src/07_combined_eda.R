suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(gprofiler2)
  source(paste0(getwd(), "/src/functions.R"))
})


c_rse <- readRDS("./data/processed/cardo_prepped_rse.rds")
s_rse <- readRDS("./data/processed/shaw_prepped_rse.rds")
t_rse <- readRDS("./data/processed/tanaka_prepped_rse.rds")

meta_comb <- t_rse@metadata %>%
  bind_rows(c_rse@metadata, s_rse@metadata) %>%
  mutate(author = ifelse(
    BioProject == "PRJNA746989", "Cardo",
    ifelse(
      BioProject == "PRJNA1093043", "Shaw", "Tanaka"
    )
  ))

write_csv(meta_comb, file = "./results/combined_metadata.csv")

# combine cardo and shaw data, subset for nonzero var genes
nonzero_genes <- intersect(rownames(c_rse), rownames(s_rse))
cts_comb <- cbind(c_rse[nonzero_genes, ], s_rse[nonzero_genes, ])
cts_filtered <- assay(cts_comb, "salmon.merged.gene_tpm")

# convert tanaka mouse IDs to human orthologs
conv <- gorth(
  rownames(
    t_rse
  ),
  source_organism = "mmusculus",
  target_organism = "hsapiens", filter_na = TRUE, mthreshold = 1
)

# subset RSEs by genes w/ mouse-human ortholog available
ortho_genes <- intersect(rownames(cts_filtered), conv$ortholog_ensg)
ortho_filtered <- cts_filtered[rownames(cts_filtered) %in% ortho_genes, ] %>%
  rownames_to_column()

# mouse genes w/ multiple orthologs for 1 human
comb_genes <- t_rse@assays@data@listData$salmon.merged.gene_tpm %>%
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
write_csv(top_cts_filt, file = "./results/combined_counts_topvar.csv")

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
summary(var_explained_df)

# select variance cutoff to filter by for retaining PC's
var_explained_cutoff <- filter(
  var_explained_df,
  var_explained > 2 / 100
)
# filter PC's for cutoff
pc_keep <- pc_df %>%
  select(seq_len(nrow(var_explained_cutoff)))

pca_meta <- meta_comb %>%
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

# color by condition
for (i in seq_along(meta_feats)) {
  rse_pca(pca_meta,
    var_explained_cutoff,
    color_by = meta_feats[i],
    shape = meta_feats[6],
    sra = "Combined Datasets",
    filepath = paste0(
      "./results/eda/combined_datasets/combined_datasets_PCA_coloredby_",
      meta_feats[i], "_WTvsHD.pdf"
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
      "./results/eda/combined_datasets/combined_datasets_PCA_coloredby_",
      meta_feats[i], "_cells.pdf"
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
      "./results/eda/combined_datasets/combined_datasets_PCA_coloredby_",
      meta_feats[i], "_author.pdf"
    )
  )
}
