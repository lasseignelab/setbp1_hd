suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
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
parser$add_argument("-m", "--meta",
                    type = "character",
                    metavar = "PATH",
                    help = "Metadata file in csv format" # nolint
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
meta <- read.csv(args$meta, row.names = 1, stringsAsFactors = TRUE)
sra <- args$sra
# meta <- read.csv("./doc/cardo_metadata.csv", row.names = 1, stringsAsFactors = TRUE)
# rse <- readRDS("./data/cardo/star_salmon/salmon.merged.gene_counts.SummarizedExperiment.rds")
# filepath <- "./results/archive/"
# var_cutoff <- 4
# sra <- "cardo"

# create output dir if needed
filepath <- paste0(args$output_dir)
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

#### Metadata ####
# dataset-specific adjustments for easier downstream analyses

if (sra == "cardo") {
  meta_filt <- meta %>%
    select(!starts_with("DATASTORE") | c(Bytes, Consent, ReleaseDate, create_date)) %>%
    # filter(rownames(.) %in% colnames(cts_filtered)) %>%
    mutate(
      # rm whitespace
      across(everything(), ~ str_replace_all(., "\\s+", "_")) ,
      cells = ifelse(cell_type == "hESC-derived_neurons", "Neuron", "NPC"),
      sample_type = paste0(cells, "_", genotype),
      sample = paste0(genotype, "_", differentiation_days),
      # broad identifier of SETBP1-KOs vs control
      condition = ifelse(
        # grepl("control", sample, ignore.case = TRUE), "Control", ifelse(
          grepl("WT", sample, ignore.case = TRUE), "Control", "SETBP1-HD"
        )) %>%
    rownames_to_column() %>%
    group_by(sample) %>%
    mutate(rep = seq_along(sample)) %>%
    column_to_rownames()
  
  print(meta_filt)
}

if (sra == "shaw") {
  meta_filt <- meta %>%
    select(!starts_with("DATASTORE") | c(Bytes, Collection_Date, Consent, ReleaseDate, create_date)) %>%
    filter(rownames(.) %in% colnames(rse)) %>%
    mutate(
      # rm whitespace
      across(everything(), ~ str_replace_all(., "\\s+", "_")) ,
      cells = ifelse(cell_type == "neural_progenitor_cells", "NPC", "IPSC"),
      sample = paste0(cells, "_", genotype),
      ACMG = ifelse(genotype == "WT", "WT", ifelse(genotype == "VUS2", "VUS", "PATH")),
      # broad identifier of SETBP1-KOs vs control
      condition = ifelse(
        # grepl("control", sample, ignore.case = TRUE), "Control", ifelse(
        grepl("WT", sample, ignore.case = TRUE), "Control", "SETBP1-HD"
      )) %>%
    rownames_to_column() %>%
    group_by(sample) %>%
    mutate(rep = seq_along(sample)) %>%
    column_to_rownames()
  
  print(meta_filt)
}

if (sra == "tanaka") {
  meta_filt <- meta %>%
    select(!starts_with("DATASTORE") | c(Bytes, Collection_Date, Consent, ReleaseDate, create_date)) %>%
    mutate(
      # rm whitespace
      across(everything(), ~ str_replace_all(., "\\s+", "_")),
      # rm prefix from samples
      sample = str_replace(`Sample.Name`, "^RNA-seq,_LSK,_", ""),
      rep = as.integer(str_extract(sample, "[0-9]+$")),
      #source_name = ifelse(`Sample.Name` starts_with("RNA"), "LSK", "LK"),
      source_name = ifelse(grepl("^RNA-seq,_LSK,_", `Sample.Name`), "LSK", "LS"),
      # split sample names and rep numbers
      sample = str_replace(sample, "-", ""),
      sample_name = str_extract(sample, ".*(?=_[0-9]+$|[0-9]+$)"),
      sample_type = ifelse(sample_name == "KO", "Setbp1_KO", ifelse(sample_name == "control", "Control", sample_name)),
      # broad identifier of SETBP1-KOs vs control
      condition = ifelse(
        grepl("control", sample, ignore.case = TRUE), "Control", "SETBP1-HD" #ifelse(
        # grepl("WT", sample, ignore.case = TRUE), "Control", "SETBP1-HD"
      )
    )
  
  print(meta_filt)
}

# save filtered/refined metadata
write_csv(meta_filt, paste0(filepath, sra, "_filtered_meta.csv"))

#### Count Data ####
##### Remove Ensembl Versions from Gene IDs #####
# ens_vers <- rownames(rse)
# ens <- str_replace_all(ens_vers, "\\..*", "")
print(paste0("Removing Ensembl versions from Gene ID's... ", head(rownames(rse))))

# update and save RSE
rownames(rse) <- gsub("\\.\\d+$", "", rownames(rse))
print(paste0("Updating rownames of RSE...", str(rse)))
rse@metadata <- meta_filt
print(paste0("Adding metadata to RSE...", str(rse)))
write_rds(rse, paste0(filepath, sra, "_prepped_rse.rds"))

##### Detect/Filter Gene Counts #####
# filtering for 0 variance genes
cts_filtered <- calc_nonzero(rse@assays@data@listData[["salmon.merged.gene_tpm"]])

write_csv(cts_filtered, paste0(filepath, sra, "_filtered_cts.csv"))