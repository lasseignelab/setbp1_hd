suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(gprofiler2)
  source(paste0(getwd(), "/src/functions.R"))
})

# t <- list.dirs("star_salmon")
# file_cts <- list.files("data/processed", pattern = "_prepped_rse.rds",
                       # full.names = TRUE)

# t <- lapply(file_cts, read_rds)

#comb_cts <- bind_cols(t @assays@data@listData[["salmon.merged.gene_tpm"]], s_rse@assays@data@listData[["salmon.merged.gene_tpm"]])
# # read in all metadata and gene counts
# c_meta <- read.csv("./doc/cardo_metadata.csv", row.names = 1, stringsAsFactors = TRUE)
# s_meta <- read.csv("./doc/shaw_metadata.csv", row.names = 1, stringsAsFactors = TRUE)
c_rse <- readRDS("./data/processed/cardo_prepped_rse.rds")
s_rse <- readRDS("./data/processed/shaw_prepped_rse.rds")
t_rse <- readRDS("./data/processed/tanaka_prepped_rse.rds")

meta_comb <- t_rse@metadata %>%
  mutate(cells = "HSC", rep = as.integer(rep)) %>%
  # mutate(rep = as.integer(rep)) %>%
  bind_rows(c_rse@metadata, s_rse@metadata)  %>%
  mutate(author = ifelse(BioProject == "PRJNA746989", "Cardo", 
                         ifelse(BioProject == "PRJNA1093043", "Shaw", "Tanaka")))

write_csv(meta_comb, file = "./results/combined_metadata.csv")

# # filepath <- "./results/archive/"
# 
# # combine tanaka data from 2 studies
# tanaka1_metadata <- read_csv("doc/tanaka1_metadata.csv")
# tanaka2_metadata <- read_csv("doc/tanaka2_metadata.csv")
# 
# tanaka_metadata <- bind_rows(tanaka1_metadata, tanaka2_metadata) %>%
#   mutate(
#     birth_date = ifelse(is.na(birth_date) == TRUE, "unknown", birth_date),
#     Collection_Date = ifelse(is.na(Collection_Date) == TRUE, "unknown", "2021-11-22")
#   )
#   #replace_na(list(birth_date = "unknown", Collection_Date = "unknown"))
#   
# write_csv(tanaka_metadata, "./doc/tanaka_metadata.csv")
# 
# combine cardo and shaw data
#cts_comb <- bind_cols(c_rse@assays@data@listData[["salmon.merged.gene_tpm"]], s_rse@assays@data@listData[["salmon.merged.gene_tpm"]])
# cts_comb <- bind_cols(assay(c_rse, "salmon.merged.gene_tpm"), assay(c_rse, "salmon.merged.gene_tpm"))
cts_comb <- cbind(c_rse, s_rse)
cts_filtered <- calc_nonzero(assay(cts_comb, "salmon.merged.gene_tpm"))

# convert tanaka mouse IDs to human orthologs
conv <- gorth(rownames(t_rse), source_organism = "mmusculus", target_organism = "hsapiens", filter_na = TRUE, mthreshold = 1)
# rse@elementMetadata@listData$human_gene_id <- conv$ortholog_ensg
# rse@elementMetadata@listData$human_gene_name <- conv$ortholog_name
# rownames(t_rse) <- conv$ortholog_ensg


# subset RSEs by genes w/ mouse-human ortholog available
#ortho_genes <- intersect(rownames(cts_comb), rownames(t_rse))
ortho_genes <- intersect(rownames(cts_filtered),conv$ortholog_ensg)
ortho_filtered <- cts_filtered[rownames(cts_filtered) %in% ortho_genes,] %>%
  rownames_to_column()

# mouse genes w/ multiple orthologs for 1 human
comb_genes <- t_rse@assays@data@listData$salmon.merged.gene_tpm %>%
  rownames_to_column() %>%
  mutate(rowname = gsub("\\.\\d+$", "", rowname)) %>%
  inner_join(conv, by = c("rowname" = "input")) %>%
  inner_join(ortho_filtered, by = c("ortholog_ensg" = "rowname"))

comb_cts <- comb_genes %>%
  select(ortholog_ensg | starts_with("SRR")) #%>%
  #column_to_rownames(var = "ortholog_ensg")
  
# human_rse <- subset(cts_comb, subset = rownames(cts_comb) %in% ortho_genes)
# mouse_rse <- subset(t_rse, subset = rownames(t_rse) %in% ortho_genes)

# full_cts <- bind_cols(assay(cts_comb, "salmon.merged.gene_tpm"), assay(mouse_rse, "salmon.merged.gene_tpm"))

# top var human genes
top_cts <- comb_cts[
  c(
    order(
      apply(comb_cts, 1, var),
      decreasing = TRUE
    )[1:5000]
  ),
]


top_cts_filt <- top_cts %>%
  # remove non-unique genes
  filter(!ortholog_ensg %in% c("ENSG00000089127", "ENSG00000099804", "ENSG00000100450", "ENSG00000110077", "ENSG00000110079", "ENSG00000116721", "ENSG00000119673", 
                              "ENSG00000136783", "ENSG00000142089", "ENSG00000157601", "ENSG00000160255", "ENSG00000167914", "ENSG00000170542", "ENSG00000172123", 
                              "ENSG00000172179", "ENSG00000198721", "ENSG00000291309")) %>%
  column_to_rownames("ortholog_ensg")

### SAVE
write_csv(top_cts_filt, file = "./results/combined_counts_topvar.csv")

# calculate PC's
pc <- prcomp(t(top_cts_filt), scale = TRUE)
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
  select(1:nrow(var_explained_cutoff))

pca_meta <- meta_comb %>%
  # filter metadata for features of interest
  select(BioProject, cells, genotype, Instrument, rep, condition, author) %>%
  # select_if(function(col) length(unique(col)) < length(col)) %>%
  bind_cols(pc_keep, .) %>%
  mutate(
    Organism = ifelse(BioProject %in% c("PRJNA898600", "PRJNA974068"), "Mouse", "Human"),
    across(where(is.character), as.factor))

meta_feats <- setdiff(colnames(pca_meta), colnames(pc_keep))

# color by condition
for (i in seq_along(meta_feats)){
  rse_pca(pca_meta, 
          var_explained_cutoff,
          color_by = meta_feats[i],
          shape_by = pca_meta$condition,
          sra = "Combined Datasets",
          filepath = paste0("./results/eda/combined_datasets/combined_datasets_PCA_coloredby_", meta_feats[i], "_WTvsHD.pdf")
  )
}
# meta_feats2 <- meta_feats[1:2]
# color by cells
for (i in seq_along(meta_feats)){
  rse_pca(pca_meta, 
          var_explained_cutoff,
          color_by = meta_feats[i],
          shape_by = pca_meta$cells,
          sra = "Combined Datasets",
          filepath = paste0("./results/eda/combined_datasets/combined_datasets_PCA_coloredby_", meta_feats[i], "_cells.pdf")
  )
}

