suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(limma)
  library(statmod)
  library(edgeR)
})

# create parser object
parser <- ArgumentParser()

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

rse <- readRDS(args$rse)
# meta <- rse@metadata
sra <- args$sra

# create output dir if needed
filepath <- paste0(args$output_dir, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

# filepath <- "./results/archive/"

# rse <- readRDS("./data/processed/cardo_prepped_rse.rds")
# sra <- "cardo"

# rse <- readRDS("./data/processed/shaw_prepped_rse.rds")
# sra <- "shaw"

# sra <- "tanaka1"
# rse <- readRDS("./data/processed/tanaka1_prepped_rse.rds")

# sra <- "tanaka2"
# rse <- readRDS("./data/processed/tanaka2_prepped_rse.rds")

cts <- "salmon.merged.gene_counts"

# create DGE list
meta <- rse@metadata %>%
  # ensure metadata samples match samples for counts
  filter(rownames(.) %in% colnames(rse)) %>%
  mutate(
    # sample = gsub("-", "_", sample_name),
    # might need to update for shaw and cardo!! double check!!
    sample = gsub("-", "_", sample),
    condition = gsub("-", "_", condition),
    across(where(is.character), as.factor)
    )

if(sra %in% c("tanaka1", "tanaka2")){
  # rename sample column
  # meta$sample <- str_to_title(meta$sample_name)
  meta$sample <- meta$condition
}
# fix for tanaka
# meta <- rse@metadata %>%
#   mutate(
#     sample = paste0(condition, sep = "_", source_name),
#     sample = gsub("-", "_", sample),
#     condition = gsub("-", "_", condition),
#     across(where(is.character), as.factor)
#   )
dge <- DGEList(counts = assay(rse, paste0(cts)), group = meta$sample, genes = as.data.frame(rse@elementMetadata@listData))

# create design matrix
# des <- model.matrix(~0 + setbp1, data = tt)
des <- model.matrix(~0 + sample, data = meta)
colnames(des) <- levels(meta$sample)
print(paste0("Number of samples per condition"))
table(meta$sample)

# filter genes with low counts
keep <- filterByExpr(dge, design = des)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# normalize and mean-variance transform with voom
dge <- calcNormFactors(dge)

if(sra == "cardo"){
  # contrast pre vs post condition
  # contr <- makeContrasts(
  #   LOWvsHI = setbp1low - setbp1high,
  #   levels = des
  # )
  contr <- makeContrasts(
    KOvsWT_15 = SETBP1_KO_15 - WT_15,
    KOvsWT_21 = SETBP1_KO_21 - WT_21,
    KOvsWT_34 = SETBP1_KO_34 - WT_34,
    KOvsWT_NPC = (SETBP1_KO_15 + SETBP1_KO_21)/2 - (WT_15 + WT_21)/2,
    KOvsWT_all = (SETBP1_KO_15 + SETBP1_KO_21 + SETBP1_KO_34)/3 - (WT_15 + WT_21 + WT_34)/3,
    levels = des
  )
}
if(sra == "shaw"){ 
  # contrast pre vs post condition
  contr <- makeContrasts(
    # individual variants vs wt
    VUSvsWT = NPC_VUS2 - NPC_WT,
    PATH2vsWT = NPC_PATH2 - NPC_WT,
    PATH3vsWT = NPC_PATH3 - NPC_WT,
    # grouped comparisons
    # NPC_PATH = (NPC_PATH2 + NPC_PATH3)/2 - (NPC_WT),
    # IPSC_PATH = (IPSC_PATH2 + IPSC_PATH3)/2 - (IPSC_WT),
    NPC_HD = (NPC_PATH2 + NPC_PATH3 + NPC_VUS2)/3 - (NPC_WT),
    IPSC_HD = (IPSC_PATH2 + IPSC_PATH3 + IPSC_VUS2)/3 - (IPSC_WT),
    # compare heterogeneity of the 2 different pathogenic variants
    PATH_vs_PATH = (NPC_PATH2 + IPSC_PATH2)/2 - (NPC_PATH3 + IPSC_PATH3)/2,
    levels = des
    )
}
if(sra %in% c("tanaka1", "tanaka2")){
  # contrast pre vs post condition
  contr <- makeContrasts(
    HDvsWT = SETBP1_HD - Control,
    # HDvsWT_LS = SETBP1_HD_LS - Control_LS,
    # HDvsWT_LSK = SETBP1_HD_LSK - Control_LSK,
    # HDvsWT_all = (SETBP1_HD_LS + SETBP1_HD_LSK)/2 - (Control_LS + Control_LSK)/2,
    levels = des
  )
}
  # limma voom
  pdf(paste0(filepath, "limmavoom_mean_var_trend.pdf"))
  v <- voom(dge, des, plot = TRUE)
  dev.off()

  # lm and empirical bayes
  fit <- lmFit(v, des)
  fit2 <- contrasts.fit(fit, contrasts = contr)
  efit <- eBayes(fit2, robust = FALSE)

  # prepare empty list for labeled DEGs
  deg_list <- vector("list", length = ncol(contr))
  names(deg_list) <- colnames(contr)

  # compute DEGs, plot, and save for each contrast
  for(i in seq_along(colnames(contr))) {
  print(paste0(colnames(contr)[i]))
  contrast <- colnames(contr)[i]
  
  res <- topTable(efit, coef = contrast, number = Inf, adjust = "BH")
  # print(paste0("Glimpse at DEG results: ", head(res)))
  
  # prep data and label for conrast and significance
  degs <- res %>%
    select(!c(transcript_id, gene_id)) %>%
    rownames_to_column(var = "gene") %>%
    mutate(contrast = paste0(contrast)) %>%
    # add significance labels
    mutate(significance = as.factor(case_when(
      logFC > 0.5 & adj.P.Val < 0.05 ~ "Up",
      logFC < -0.5 & adj.P.Val < 0.05 ~ "Down",
      # logFC > 0.5 & adj.P.Val < 0.1 ~ "Up",
      # logFC < -0.5 & adj.P.Val < 0.1 ~ "Down",
       .default = "NonSig"
    )))
  print(paste0("DEGs summary: ")) 
  print(table(degs$significance))
  
  # add degs to list
  deg_list[[i]] <- degs
  
  # save results
  write_csv(res, file = paste0(filepath, contrast, "_unfiltered_degs.csv"))
}

# combine degs across contrasts
degs_comb <- bind_rows(deg_list)
sig_degs_comb <- degs_comb %>%
  filter(significance != "NonSig") %>%
  mutate(se_logFC = abs(logFC / t))

# save labeled and sig-filtered results
write_csv(degs_comb, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
write_csv(sig_degs_comb, file = paste0(filepath, "combined_contrasts_sig_degs.csv"))

# volcano plots across contrasts
vol_plot <- degs_comb %>%
  # filter out highest p-adj genes for easier plotting
  #filter(adj.P.Val < 0.5) %>%
  ggplot(aes(y = -log10(`adj.P.Val`), x = logFC, color = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +
  facet_wrap(~ contrast)

vol_plot
ggsave(filename = paste0(filepath,  "degs_volcanoplot.png"))

# DEG counts per contrast
degs_comb %>%
  ggplot(aes(y = contrast, fill = significance)) +
  geom_bar(position = position_dodge()) +
  scale_x_log10()
ggsave(filename = paste0(filepath,  "degs_summary_bar.png"))

# TPM of top DEGs
top_degs <- sig_degs_comb %>%
  slice_max(n = 50, order_by = abs(logFC))

top_cts <- assay(rse, "salmon.merged.gene_tpm") %>%
  filter(rownames(.) %in% top_degs$gene) %>%
  mutate_if(is.numeric, ~ log2(. + 1)) %>%
  rownames_to_column(var = "gene") %>%
  left_join(top_degs[,c("gene", "gene_name")], by = "gene") %>%
  select(!gene) %>%
  distinct()
  # column_to_rownames(var = "gene_name")
  # bind_cols(top_degs[,c("gene", "gene_name")], by)
top_mat <- as.matrix(top_cts[, !names(top_cts) %in% c("gene_name")]) 

# Set column names to sample names
colnames(top_cts) <- c(paste0(meta$sample), "gene_name")

# Create heatmap
png(paste0(filepath, "top_degs_heatmap.png"), width = 22, height = 24, units = "cm", res = 300)
set.seed(42)
sample_anno <- HeatmapAnnotation(df = meta$sample)

Heatmap(top_mat,# ,colnames(top_cts) %in% c("gene_name"))],
        name = "log2(TPM)",
        col = colorRamp2(c(min(top_mat), 0, max(top_mat)), 
                         c("blue", "white", "red")),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
  row_labels = top_cts$gene_name,
  bottom_annotation = sample_anno,
        # show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10),
        # row_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "Expression", 
                                    title_gp = gpar(fontsize = 10, fontface = "bold"),
                                    labels_gp = gpar(fontsize = 8)))

dev.off()



# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
