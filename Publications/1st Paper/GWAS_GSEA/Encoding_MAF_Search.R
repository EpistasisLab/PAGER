### R
### PJF
### 2/23/24

### This script will match the genotype scores from Additive to PAGER and EDGE for putative QTL (or any list of SNPs) to determine the encoding scores to detect patterns for beta coefficient correction and inheritance model conformity. It also calculates the MAF for putative QTL

library(data.table)
library(tidyverse)

setwd("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/GWAS")

# Load Genotype Files
Additive <- fread("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/RatData-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/Geno_noPAGER_nophen_no_extra_columns.csv")
PAGER <- fread("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/RatData-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/Geno_PAGER_nophen_no_extra_columns.csv")
EDGE <- fread("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/RatData-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/OtherPhen_PAGER_EDGE_files/Geno_EDGE_BW_Encoded_with_NA.csv")

# Check PAGER encodings of all unique putative QTL SNP locations
# Load in SNP list
snpList <- c("chr1.106956497", "chr1.160493012", "chr1.160493164", "chr1.185170500", "chr1.185730317", "chr1.281420356", "chr1.281489331", "chr1.281509176", "chr1.282070632", "chr2.65816485", "chr2.69780151", "chr2.72967800", "chr3.95378412", "chr3.136021511", "chr3.136707086", "chr3.137286589", "chr3.137335822", "chr5.50933779", "chr6.28560776", "chr7.24971798", "chr7.25030336", "chr7.34874677", "chr7.34917462", "chr8.81316464", "chr10.84021443", "chr10.110912475", "chr12.5046459", "chr12.5747404", "chr12.5782645", "chr13.53419980", "chr18.27223157", "chr18.27355039") 

# Loop through SNPs and retrieve PAGER values
for(snp in snpList) {
  if(snp %in% colnames(Additive)) {
    # For each genotype (0, 1, 2), find the first occurrence's row index in Additive
    genotypeIndexes <- map_dbl(c(0, 1, 2), ~which(Additive[[snp]] == .x)[1])
    
    # Map these indexes to PAGER encodings
    pagerEncodings <- map_chr(genotypeIndexes, ~as.character(PAGER[[snp]][.x]))
    
    # Print the mappings for this SNP
    cat(sprintf("For SNP %s: 0 = '%s', 1 = '%s', and 2 = '%s'\n", 
                snp, pagerEncodings[1], pagerEncodings[2], pagerEncodings[3]))
  } else {
    cat(sprintf("SNP %s not found in the Additive dataframe.\n", snp))
  }
}

# Loop through SNPs and retrieve EDGE values for heterozygotes
for(snp in snpList) {
  if(snp %in% colnames(Additive) && snp %in% colnames(EDGE)) {
    # Find the first occurrence's row index of heterozygote (1) in Additive
    heterozygoteIndex <- which(Additive[[snp]] == 1)[1]
    
    # Check if a heterozygote index was found
    if(!is.na(heterozygoteIndex)) {
      # Retrieve the EDGE encoding for the heterozygote
      edgeEncoding <- EDGE[[snp]][heterozygoteIndex]
      
      # Print the mapping for this SNP heterozygote
      cat(sprintf("For SNP %s, heterozygote (1) = '%s'\n", snp, edgeEncoding))
    } else {
      # If no heterozygote is found in Additive for this SNP
      cat(sprintf("No heterozygote (1) found for SNP %s in Additive.\n", snp))
    }
  } else {
    if(!(snp %in% colnames(Additive))) {
      cat(sprintf("SNP %s not found in the Additive dataframe.\n", snp))
    } else {
      cat(sprintf("SNP %s not found in the EDGE dataframe.\n", snp))
    }
  }
}

# Calculate minor allele frequencies of all SNP QTL
mafList <- sapply(snpList, function(snp) {
  if(snp %in% colnames(Additive)) {
    genotypes <- Additive[[snp]] # Access the SNP column
    nonNA_genotypes <- genotypes[!is.na(genotypes)] # Exclude NA values
    if(length(nonNA_genotypes) == 0) return(NA) # Return NA if all values are NA
    totalAlleles <- 2 * length(nonNA_genotypes) # Adjusted for non-NA genotypes
    totalAltAlleles <- 2 * sum(nonNA_genotypes == 2, na.rm = TRUE) + sum(nonNA_genotypes == 1, na.rm = TRUE)
    freqAlt <- totalAltAlleles / totalAlleles
    maf <- min(freqAlt, 1 - freqAlt)
    return(maf)
  } else {
    return(NA) # Return NA if SNP not found
  }
}, USE.NAMES = TRUE)

# Print the MAF for each SNP
print(mafList)
