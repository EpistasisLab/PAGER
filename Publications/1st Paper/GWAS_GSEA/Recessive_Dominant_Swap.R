library(data.table)

setwd("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK")

Dominant_GWAS <- fread("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/Geno_noPAGER_nophen_no_extra_columns.csv")

# Specify the columns you want to check and replace
columns_to_check <- 1:128401

# Check if all three values 0, 1, and 2 are present in the specified columns
columns_to_replace <- columns_to_check[sapply(Dominant_GWAS[, ..columns_to_check, with = FALSE], function(x) 0 %in% x & 1 %in% x & 2 %in% x)]

## Create Recessive dataset
# Replace 0 with 1 in the selected columns
Dominant_GWAS[, (columns_to_replace) := lapply(.SD, function(x) ifelse(x == 0, 1, x)), .SDcols = columns_to_replace]

# Replace 2 with 0 in the selected columns
Dominant_GWAS[, (columns_to_replace) := lapply(.SD, function(x) ifelse(x == 2, 0, x)), .SDcols = columns_to_replace]

fwrite(Dominant_GWAS, "Geno_Dominant_nophen_no_extra_columns.csv", row.names = FALSE)

## Create Dominant dataset
Recessive_GWAS <- fread("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/Rat Data-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/Geno_noPAGER_nophen_no_extra_columns.csv")

# Convert Additive to Dominant encoding
Recessive_GWAS[, (columns_to_replace) := lapply(.SD, function(x) ifelse(x == 2, 1, x)), .SDcols = columns_to_replace]

fwrite(Recessive_GWAS, "Geno_Recessive_nophen_no_extra_columns.csv", row.names = FALSE)
