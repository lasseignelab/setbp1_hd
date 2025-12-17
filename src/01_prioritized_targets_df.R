suppressPackageStartupMessages({
 library(readr)
 library(dplyr)
 library(readxl)
})

# BrainSpan data 
rows_meta <- read_csv()
cols_meta <- read_csv()
raw_expression <- read_csv()

# SETBP1 prioritized targets list 
prioritized_targets_df <-read_csv()
prioritized_targets <- unique(prioritized_targets_df[[1]])

# Defining 24-37 pcw into years
24_weeks <- (24-40)/52
37_weeks <- (37-40)/52

24_to_37_weeks <-cols_meta$age_years >= 24_weeks & cols_meta$age_years <= 37_weeks

cat("Number of samples in 24-37 gestation window:", sum(24_to_37_weeks), "\n")
if (sum(24_to_37_weeks) < 3) stop ("Not enough samples")

age_window <- cols_meta$age_years[24_to_37_weeks]




