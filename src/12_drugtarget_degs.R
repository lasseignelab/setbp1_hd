# script for analyzing DEGs of drug targets from signature reversion results
set.seed(42)

suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(UpSetR)
  library(ComplexHeatmap)
  library(statmod)
  library(gprofiler2)
  library(circlize)
  library(tidyverse)
  library(readxl)
})
###### SET UP ######
# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--indir",
  type = "character",
  help = "Path to the input data directory to access validation data"
)
parser$add_argument("-d", "--degs",
  type = "character",
  help = "Parent directory containing each dataset's DEG CSV files"
)
parser$add_argument("-o", "--outdir",
  type = "character",
  help = "Path to the directory to read in signature results and write target DEG results" # nolint
)
# get command line options, if help option encountered print help and exit
args <- parser$parse_args()
indir <- args$indir
degs <- args$degs
outdir <- args$outdir

# create output dir if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# helper function to read in DEGs
read_deg <- function(dataset_name) {
  file <- file.path(degs, dataset_name, "significant_degs.csv")
  if (!file.exists(file)) stop("Missing file: ", file)
  read_csv(file) %>% mutate(dataset = dataset_name)
}

# start timer
ptm <- proc.time()

###### READ IN DATA ######
drugs_brain <- read_rds(
  paste0(outdir, "lincs_anno_braincell_results.rds")
)
drugs_setbp1 <- read_csv(
  paste0(outdir, "setbp1_drug_targets.csv")
)

# Wong et al data
wong_fibro <- read_excel(paste0(
  indir, "wong/SourceData/SourceData6_forFIg5c_DEG.xlsx"
), sheet = 3, skip = 1) %>%
  select(
    gene = ensemblID, logFC = log2FoldChange,
    gene_name = symbol, adj.P.Val = padj
  ) %>%
  mutate(contrast = "fibro", dataset = "wong")

wong_neu <- read_excel(paste0(
  indir, "wong/SourceData/SourceData10_forFig7b_DEG.xlsx"
), sheet = 3, skip = 1) %>%
  select(
    gene = `...1`, logFC = log2FoldChange,
    gene_name = symbol, adj.P.Val = padj
  ) %>%
  mutate(contrast = "neu", dataset = "wong")

# combine validation DEGs and cardo/shaw DEGs
all_degs <- bind_rows(
  read_deg("cardo"),
  read_deg("shaw"),
  read_deg("tanaka1"),
  read_deg("tanaka2"),
  wong_neu,
  wong_fibro
)

### SIMPLIFIED DRUG TABLE ###
brain_drugs_filt <- drugs_brain %>%
  dplyr::select(
    pert, pref_name, max_phase, cell,
    WTCS, WTCS_Pval, mergeTargets, Target_pathway
  ) %>%
  filter(!is.na(pref_name))
write_csv(brain_drugs_filt, paste0(outdir, "brain_drugs_filtered.csv"))

###### FIND DRUG TARGET DEGS ######
drugs_brain_targets <- drugs_brain %>%
  select(
    pert, pert_iname, pref_name, WTCS, WTCS_Pval, mergeTargets
  ) %>%
  separate_longer_delim(mergeTargets, delim = "; ") %>%
  filter(!is.na(mergeTargets))

# exclude redundant/biased contrasts
contrasts_skip <- c(
  "KOvsWT_NPC", "KOvsWT_all", "NPC_HD", "PATH_vs_PATH"
)

target_degs <- all_degs %>%
  filter(gene_name %in% drugs_brain_targets$mergeTargets)
message(
  "Number of unique target DEGs: ", length(unique(target_degs$gene_name))
)
message(
  "Target DEGs that are also SETBP1 targets: ",
  intersect(target_degs$gene_name, drugs_setbp1$mergeTargets)
)
message(
  "SETBP1 targets that are not DEGs: ",
  setdiff(drugs_setbp1$mergeTargets, target_degs$gene_name)
)

# save consensus DEGs
write_csv(target_degs, paste0(outdir, "drug_target_degs.csv"))

target_degs %>%
  filter(gene_name %in% drugs_setbp1$mergeTargets) %>%
  ggplot(aes(
    x = -log10(adj.P.Val),
    y = logFC,
    colour = contrast, shape = gene_name
  )) +
  geom_point()
ggsave(paste0(outdir, "target_exp_contrasts.png"))

target_degs_wide <- target_degs %>%
  unite("contrast_dataset", contrast, dataset) %>%
  select(gene, gene_name, logFC, contrast_dataset) %>%
  distinct() %>%
  pivot_wider(
    names_from = contrast_dataset, values_from = logFC, values_fill = 0
  ) %>%
  # rm characters "NA" -- not in NA class
  filter(gene_name != "NA")

str(target_degs_wide)

### prepare Heatmap parameters ###
message("Preparing heatmap parameters...")
# color scale for logFC
max_lfc <- max(abs(target_degs$logFC))
col_fun <- circlize::colorRamp2(
  c(-max_lfc, 0, max_lfc), c("blue", "white", "red")
)

# column annotations
mat_cols <- colnames(target_degs_wide)[3:ncol(target_degs_wide)]

# Define groups using case_when logic
col_groups <- case_when(
  grepl("_cardo$", mat_cols) ~ "cardo",
  grepl("_shaw$", mat_cols) ~ "shaw",
  TRUE ~ "validation"
)
col_groups

# Define specific colors
col_anno_colors <- c(
  "cardo"      = "#4DAF4A", # Black
  "shaw"       = "#377EB8", # Blue
  "validation" = "#984EA3" # Purple
)

column_ha <- HeatmapAnnotation(
  Source = col_groups,
  col = list(Source = col_anno_colors),
  show_annotation_name = FALSE
)

#--- Row Annotations ---#
row_groups <- target_degs_wide %>%
  select(gene, gene_name) %>%
  mutate(
    SETBP1_target = ifelse(
      gene_name %in% unique(drugs_setbp1$mergeTargets), TRUE, FALSE
    )
  ) %>%
  column_to_rownames(var = "gene")

gene_ha <- rowAnnotation(
  SETBP1_target = row_groups$SETBP1_target,
  col = list(SETBP1_target = c("TRUE" = "steelblue", "FALSE" = "white"))
)

### plot heatmap ###
message("Plotting all target DEGs heatmap...")

heatmap_obj <- Heatmap(
  as.matrix(target_degs_wide[, 3:length(target_degs_wide)]),
  name = "Log2(Fold Change)",
  col = col_fun,
  row_labels = target_degs_wide$gene_name,
  row_names_gp = gpar(fontsize = 4.6),
  bottom_annotation = column_ha,
  right_annotation = gene_ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  border = TRUE,
  column_title = "Significant Drug Target DEGs Across Studies"
)

# save
png(paste0(
  outdir, "target_exp_heatmap.png"
), width = 20, height = 55, units = "cm", res = 300)
draw(heatmap_obj)
dev.off()

### DRUG TARGETS FILTERED FOR APPROVED DRUGS ###
approved_drugs <- drugs_brain %>%
  filter(max_phase == 4) %>%
  separate_longer_delim(mergeTargets, delim = "; ")
message("Approved drug drug targets: ", unique(approved_drugs$mergeTargets))

target_degs_approved <- target_degs_wide %>%
  filter(gene_name %in% approved_drugs$mergeTargets)

#--- Row Annotations ---#
row_groups <- target_degs_approved %>%
  select(gene, gene_name) %>%
  mutate(SETBP1_target = ifelse(
    gene_name %in% unique(drugs_setbp1$mergeTargets), TRUE, FALSE
  )) %>%
  column_to_rownames(var = "gene")

gene_ha <- rowAnnotation(
  SETBP1_target = row_groups$SETBP1_target,
  col = list(SETBP1_target = c("TRUE" = "steelblue", "FALSE" = "white"))
)
### plot heatmap ###
message("Plotting approved drugs target DEGs heatmap...")

heatmap_obj <- Heatmap(
  as.matrix(target_degs_approved[, 3:length(target_degs_approved)]),
  name = "Log2(Fold Change)",
  col = col_fun,
  row_labels = target_degs_approved$gene_name,
  row_names_gp = gpar(fontsize = 4.6),
  bottom_annotation = column_ha,
  right_annotation = gene_ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  border = TRUE,
  column_title = "Significant Approved Drug Target DEGs Across Studies"
)

# save
png(paste0(
  outdir, "approveddrugs_target_exp_heatmap.png"
), width = 20, height = 55, units = "cm", res = 300)
draw(heatmap_obj) # , column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()


#### END ####

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
