## 1/4/24
## PJF
## R
## Determining the 5% percentile of permutated Pvals for each encoding using an extreme value distribution

setwd("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/lowest_pval_results")

library(evd)

# Import data
additive <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/lowest_pval_results/additive_gwas/bmiwtail/lowest_p_score_values_additive_gwas.csv")

dominant <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/lowest_pval_results/dominant_additive_gwas/bmiwtail/lowest_p_score_values_dominant_with_additive_GRMs.csv")

recessive <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/lowest_pval_results/recessive_additive_gwas/lowest_p_score_values_recessive_with_additive_GRMs.csv")

EDGE <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/lowest_pval_results/edge_gwas/lowest_p_score_values_edge_gwas_second_permutations.csv")

PAGER <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/lowest_pval_results/pager_gwas/lowest_p_score_values_pager_gwas_second_permutations.csv")

## Additive
# Assuming you have a vector of your lowest p-values from each permutation
additive_lowp <- additive$Lowest_p_score  # Replace with your data

# Histogram of your data
hist(additive_lowp, breaks = 30, probability = TRUE, col = "lightblue", main = "Histogram with Fitted Gumbel Distribution", xlab = "P-Values")

# Fit the data to a Gumbel distribution
additive_fit <- fgev(additive_lowp, std.err = FALSE)

# Check the summary of the fit
summary(additive_fit)
additive_fit$estimate["shape"] # 1.396947e-105

# To find the significance threshold, you would use the qgev function
desired_percentile <- 0.95  # For the 95th percentile
add_significance_threshold <- qgev(1 - desired_percentile, loc = additive_fit$estimate["loc"], scale = additive_fit$estimate["scale"], shape = additive_fit$estimate["shape"])
-log10(add_significance_threshold) # 5.45

## Dominant
# Assuming you have a vector of your lowest p-values from each permutation
dominant_lowp <- dominant$Lowest_p_score  # Replace with your data

# Histogram of your data
hist(dominant_lowp, breaks = 30, probability = TRUE, col = "lightblue", main = "Histogram with Fitted Gumbel Distribution", xlab = "P-Values")

# Fit the data to a Gumbel distribution
dominant_fit <- fgev(dominant_lowp, std.err = FALSE)

# Check the summary of the fit
summary(dominant_fit)
dominant_fit$estimate["shape"] # 2.628208e-92

# To find the significance threshold, you would use the qgev function
desired_percentile <- 0.95  # For the 95th percentile
dom_significance_threshold <- qgev(1 - desired_percentile, loc = dominant_fit$estimate["loc"], scale = dominant_fit$estimate["scale"], shape = dominant_fit$estimate["shape"])
-log10(dom_significance_threshold) # 5.35

## Recessive
# Assuming you have a vector of your lowest p-values from each permutation
recessive_lowp <- recessive$Lowest_p_score  # Replace with your data

# Histogram of your data
hist(recessive_lowp, breaks = 30, probability = TRUE, col = "lightblue", main = "Histogram with Fitted Gumbel Distribution", xlab = "P-Values")

# Fit the data to a Gumbel distribution
recessive_fit <- fgev(recessive_lowp, std.err = FALSE)

# Check the summary of the fit
summary(recessive_fit)
recessive_fit$estimate["shape"] # 2.479619e-80

# To find the significance threshold, you would use the qgev function
desired_percentile <- 0.95  # For the 95th percentile
rec_significance_threshold <- qgev(1 - desired_percentile, loc = recessive_fit$estimate["loc"], scale = recessive_fit$estimate["scale"], shape = recessive_fit$estimate["shape"])
-log10(rec_significance_threshold) # 5.44

## EDGE
# Assuming you have a vector of your lowest p-values from each permutation
EDGE_lowp <- EDGE$Lowest_p_score  # Replace with your data

# Histogram of your data
hist(EDGE_lowp, breaks = 30, probability = TRUE, col = "lightblue", main = "Histogram with Fitted Gumbel Distribution", xlab = "P-Values")

# Fit the data to a Gumbel distribution
EDGE_fit <- fgev(EDGE_lowp, std.err = FALSE)

# Check the summary of the fit
summary(EDGE_fit)
EDGE_fit$estimate["shape"] # 1.521784e-77 

# To find the significance threshold, you would use the qgev function
desired_percentile <- 0.95  # For the 95th percentile
EDGE_significance_threshold <- qgev(1 - desired_percentile, loc = EDGE_fit$estimate["loc"], scale = EDGE_fit$estimate["scale"], shape = EDGE_fit$estimate["shape"])
-log10(EDGE_significance_threshold) # 5.37

## PAGER # Setting loc manually
# Assuming you have a vector of your lowest p-values from each permutation
PAGER_lowp <- PAGER$Lowest_p_score  # Replace with your data

# Histogram of your data
hist(PAGER_lowp, breaks = 30, probability = TRUE, col = "lightblue", main = "Histogram with Fitted Gumbel Distribution", xlab = "P-Values")

# Fit the data to a Gumbel distribution
PAGER_fit <- fgev(PAGER_lowp, std.err = FALSE, loc = median(PAGER_lowp))

# Check the summary of the fit
summary(PAGER_fit)
PAGER_fit$estimate["shape"] # 1.453163e-17 

# To find the significance threshold, you would use the qgev function
desired_percentile <- 0.95  # For the 95th percentile
PAGER_significance_threshold <- qgev(1 - desired_percentile, loc = median(PAGER_lowp), scale = PAGER_fit$estimate["scale"], shape = PAGER_fit$estimate["shape"]) 
-log10(PAGER_significance_threshold) # 6.06
