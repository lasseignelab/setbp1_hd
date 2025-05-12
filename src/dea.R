suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(limma)
  library(gprofiler2)
  # source(paste0(getwd(), "/functions.R"))
  source(paste0(getwd(), "/src/functions.R"))
})


# rse <- readRDS(args$rse)
# meta <- rse@metadata
# var_cutoff <- args$var_cutoff
# sra <- args$sra

# rse <- readRDS("./data/processed/cardo_prepped_rse.rds")
filepath <- "./results/archive/"
rse <- readRDS("./data/processed/shaw_prepped_rse.rds")
cts <- "salmon.merged.gene_tpm"

meta <- rse@metadata %>%
  mutate(
    sample = gsub("-", "_", sample),
    condition = gsub("-", "_", condition),
    across(where(is.character), as.factor)
    )

condition <- meta$condition
# days <- meta$differentiation_days
baseline <- "Control"


# gene cts filtered for nonzero var
cts_filtered <- calc_nonzero(assay(rse, "salmon.merged.gene_tpm"))


# print(paste0("Ensure non 'baseline' options are all referring to post-condition (or else contrasts will be incorrect. Unique condition status options: ", unique(rse@colData@listData$sra_attribute.condition)))

  # make design
  # group <- paste0(meta$condition, meta$differentiation_days)
  # des <- model.matrix(~condition + days)
  # des <- model.matrix(~0 + condition)
  # des <- model.matrix(~0 + sample + days, data = meta)
des <- model.matrix(~0 + sample, data = meta)
  
  # contrast pre vs post condition
  contr <- makeContrasts(
    # WTvsHD = (sampleWT_15 + sampleWT_21 + sampleWT_34)/3-(sampleSETBP1_KO_15 + sampleSETBP1_KO_21 + sampleSETBP1_KO_34)/3,
    # WTvsHD = (sampleWT_15 + sampleWT_21 + sampleWT_34)-(sampleSETBP1_KO_15 + sampleSETBP1_KO_21 + sampleSETBP1_KO_34),
    # WTvsHD = conditionControl - conditionSETBP1_HD,
    # WTvsHD34 = sampleWT_34 - sampleSETBP1_KO_34,
    NPC_PATH = (sampleNPC_PATH2 + sampleNPC_PATH3)/2 - (sampleNPC_WT),
    IPSC_PATH = (sampleIPSC_PATH2 + sampleIPSC_PATH3)/2 - (sampleIPSC_WT),
    NPC_HD = (sampleNPC_PATH2 + sampleNPC_PATH3 + sampleNPC_VUS2)/3 - (sampleNPC_WT),
    IPSC_HD = (sampleIPSC_PATH2 + sampleIPSC_PATH3 + sampleIPSC_VUS2)/3 - (sampleIPSC_WT),
    levels = des
    )
  # contr <- makeContrasts(preVSpost = conditionControl-conditionSETBP1_HD, levels = des) #colnames(des))
  # contr <- makeContrasts(contrasts = c(conditionControl:days15-conditionSETBP1_HD:days15,
  #                                      conditionControl:days21-conditionSETBP1_HD:days21,
  #                                      conditionControl:days34-conditionSETBP1_HD:days34),
  #                        levels = des)
                                       #))# ctl_vs_hd = conditionControl-conditionSETBP1_HD, levels = des)
  
  # LIMMA and empirical bayes, take top 100 DEGs
  fit <- lmFit(cts_filtered, des)
  fit <- contrasts.fit(fit, contrasts = contr)
  efit <- eBayes(fit)
  # plotSA(efit)
  degs <- topTable(efit, n = 10000,  p.value = 0.5)
  
  # ens to symbol conversion
  t <- gconvert(rownames(degs))
  degs_long <- degs %>%
    rownames_to_column(var = "input") %>%
    inner_join(t, c = "input") %>%
    pivot_longer(cols = c(NPC_PATH, IPSC_PATH, NPC_HD, IPSC_HD), 
                 names_to = "Comparison", values_to = "logFC")

vol_plot <- degs_long %>%
  ggplot(aes(y = -log10(`adj.P.Val`), x = log2(logFC), color = Comparison)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed")

vol_plot
ggsave(filename = "./results/archive/DEA_shaw_volcanoplot.png")

sig_degs <- decideTests(fit)

write_rds(sig_degs, file = "./results/archive/DEA_shaw_sig_degs.rds")

common_path_degs <- which(sig_degs[,1]!=0 & sig_degs[,2])


