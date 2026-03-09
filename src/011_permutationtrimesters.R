#### Permutation Analysis of SETBP1 Targets by Trimester ####
# GOAl: Getting a closer view to determine when this change in gene expression occur. 
## NOTE: 1st = 1-14 wks pcw; 2nd = 14-27 pcw; 3rd = 28-40 pcs

# RESULTS: Largest change is shown between 1st and 2nd Semester (p < 0.05)

# Creating Sub-windows by trimester
## Grouping BrainSpan pcw into trimesters 
cols_meta <- cols_meta %>%
  mutate(
    trimester = case_when(
      age_years < (14 - 40) / 52 ~ "First_Trimester",
      age_years >= (14 - 40) / 52 & age_years < (28 - 40) / 52 ~ "Second_Trimester",
      age_years >= (28 - 40) / 52 & age_years < 0 ~ "Third_Trimester",
      TRUE ~ "Not_Prenatal"
    )
  )

# Categorizing SETBP1 target expression by trimester
tri1 <- cols_meta$trimester == "First_Trimester"
tri2 <- cols_meta$trimester == "Second_Trimester"
tri3 <- cols_meta$trimester == "Third_Trimester"

# Calculate means for SETBP1 targets 
tri1_setbp1 <- rowMeans(expr_mat[setbp1_targets_use, tri1, drop = FALSE], na.rm = TRUE)
tri2_setbp1 <- rowMeans(expr_mat[setbp1_targets_use, tri2, drop = FALSE], na.rm = TRUE)
tri3_setbp1 <- rowMeans(expr_mat[setbp1_targets_use, tri3, drop = FALSE], na.rm = TRUE)

cat("SETBP1 Target Averages by Trimester\n")
cat("First Trimester (<14 pcw):  ", mean(tri1_setbp1, na.rm = TRUE), "\n")
cat("Second Trimester (14-27 pcw):", mean(tri2_setbp1, na.rm = TRUE), "\n")
cat("Third Trimester (28-40 pcw): ", mean(tri3_setbp1, na.rm = TRUE), "\n")

# Confirming differences in averages are significant
## First vs Second
first_vs_second <- wilcox.test(tri1_setbp1, tri2_setbp1, paired = TRUE, alternative = "greater")
cat("First vs Second Trimester\n")
cat("p-value =", first_vs_second$p.value, "\n")
if (first_vs_second$p.value < 0.05) {
  cat ("First vs Second = SIGNIFICANT\n")
}

## First vs Third
first_vs_third <-  wilcox.test(tri1_setbp1, tri3_setbp1, paired = TRUE, alternative = "greater")
cat("First vs Third Trimester\n")
cat("p-value =", first_vs_third$p.value, "\n")
if (first_vs_second$p.value < 0.05) {
  cat ("First vs Third = SIGNIFICANT\n")
}

## Second vs Third
second_vs_third <- wilcox.test(tri2_setbp1, tri3_setbp1, paired = TRUE, alternative = "greater")
cat("Second vs Third Trimester\n")
cat("p-value =", second_vs_third$p.value, "\n")
if (second_vs_third$p.value < 0.05) {
  cat ("Second vs Third = SIGNIFICANT\n")
}


