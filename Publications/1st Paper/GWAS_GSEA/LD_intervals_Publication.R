### PJF
### R
### 2/21/24
### This script loops through a list of SNPs and caculates the LD interval around it using Plink.

# Load necessary library
library(data.table)

# Set up directories and project name
setwd("/Users/fredap/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/RatData-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/LD_CSA")

exp_dir <- getwd() # Current working directory

## Additive
# Read in the list of top SNPs from a CSV file
topsnps_add <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/RatData-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/LD_CSA/additive_QTLs_LD.csv")

topsnps_add$Additive <- gsub("\\.", ":", topsnps_add$Additive) # Replace "." with ":" in the first column

colnames(topsnps_add) <- c("topsnps", "trait") # Change column names to match script below

# Extract chromosome and position from the 'topsnp' column
topsnps_add$chr <- sapply(strsplit(topsnps_add$topsnp, split=":"), `[`, 1)
topsnps_add$pos <- sapply(strsplit(topsnps_add$topsnp, split=":"), `[`, 2)

# Define a function to calculate LD interval start and stop positions
LD_start_stop <- function(x, trait, chr, pos) {
  data <- fread(x, header=T, stringsAsFactors =F, select=c(3,6,7))
  setnames(data, c("snp1", "snp2", "rsquare"))
  data$dprime <- rep(NA, nrow(data))
  data <- data[which(data$rsquare >= 0.6),]
  
  # Split 'snp2' to extract chromosome and position, convert position to numeric
  data[, c("chr", "pos") := tstrsplit(snp2, ":", fixed=TRUE)]
  data$pos <- as.numeric(data$pos)
  data <- data[order(pos),]
  
  # Determine the start and stop positions of the LD interval
  LD_interval_start <- data$pos[1]
  LD_interval_stop <- data$pos[nrow(data)]
  output <- data.frame(LD_interval_start, LD_interval_stop)
  return(output)
}

# Initialize an empty dataframe to store results
out <- data.frame()

# Loop over each SNP to calculate LD intervals
for(i in 1:nrow(topsnps_add)) {
  trait <- topsnps_add$trait[i]
  topsnp <- topsnps_add$topsnp[i]
  chr <- topsnps_add$chr[i]
  pos <- topsnps_add$pos[i]
  
  # Run PLINK to compute LD, outputting results to a temporary file
  system(paste0("./plink --bfile ",exp_dir,"/P50_round2_LD_pruned_3473 --chr ",gsub("chr", "", chr)," --nonfounders --r2  --ld-snp ",topsnp," --ld-window 100000 --ld-window-kb 6000 --ld-window-r2 0 --out ",exp_dir,"/LZ_temp"),ignore.stdout = T,ignore.stderr = T,wait = T)
  
  # Calculate LD interval start and stop positions
  test <- LD_start_stop(paste0(exp_dir, "/LZ_temp.ld"), trait = trait, chr = chr, pos = pos)
  out <- rbind(out, test) # Append results to the output dataframe
}

# Bring LD interval information into the dataframe
topsnps_add$LD_start <- out$LD_interval_start
topsnps_add$LD_stop <- out$LD_interval_stop

# Write updated file to disk
write.csv(topsnps_add, "additive_QTLs_LD.csv", row.names = FALSE)

## EDGE
# Read in the list of top SNPs from a CSV file
topsnps_EDGE <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/RatData-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/LD_CSA/EDGE_QTLs_LD.csv")

topsnps_EDGE$EDGE <- gsub("\\.", ":", topsnps_EDGE$EDGE) # Replace "." with ":" in the first column

colnames(topsnps_EDGE) <- c("topsnps", "trait") # Change column names to match script below

# Extract chromosome and position from the 'topsnp' column
topsnps_EDGE$chr <- sapply(strsplit(topsnps_EDGE$topsnp, split=":"), `[`, 1)
topsnps_EDGE$pos <- sapply(strsplit(topsnps_EDGE$topsnp, split=":"), `[`, 2)

# Initialize an empty dataframe to store results
out <- data.frame()

# Loop over each SNP to calculate LD intervals
for(i in 1:nrow(topsnps_EDGE)) {
  trait <- topsnps_EDGE$trait[i]
  topsnp <- topsnps_EDGE$topsnp[i]
  chr <- topsnps_EDGE$chr[i]
  pos <- topsnps_EDGE$pos[i]
  
  # Run PLINK to compute LD, outputting results to a temporary file
  system(paste0("./plink --bfile ",exp_dir,"/P50_round2_LD_pruned_3473 --chr ",gsub("chr", "", chr)," --nonfounders --r2  --ld-snp ",topsnp," --ld-window 100000 --ld-window-kb 6000 --ld-window-r2 0 --out ",exp_dir,"/LZ_temp"),ignore.stdout = T,ignore.stderr = T,wait = T)
  
  # Calculate LD interval start and stop positions
  test <- LD_start_stop(paste0(exp_dir, "/LZ_temp.ld"), trait = trait, chr = chr, pos = pos)
  out <- rbind(out, test) # Append results to the output dataframe
}

# Bring LD interval information into the dataframe
topsnps_EDGE$LD_start <- out$LD_interval_start
topsnps_EDGE$LD_stop <- out$LD_interval_stop

# Write updated file to disk
write.csv(topsnps_EDGE, "EDGE_QTLs_LD.csv", row.names = FALSE)

## PAGER
# Read in the list of top SNPs from a CSV file
topsnps_PAGER <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/RatData-PalmerLab/Genotypes/All_GWAS_files/LDPruned_PAGER/Genotypes/LD_pruned_PLINK/QTLPrune_Sig/LD_CSA/PAGER_QTLs_LD.csv")

topsnps_PAGER$PAGER <- gsub("\\.", ":", topsnps_PAGER$PAGER) # Replace "." with ":" in the first column

colnames(topsnps_PAGER) <- c("topsnps", "trait") # Change column names to match script below

# Extract chromosome and position from the 'topsnp' column
topsnps_PAGER$chr <- sapply(strsplit(topsnps_PAGER$topsnp, split=":"), `[`, 1)
topsnps_PAGER$pos <- sapply(strsplit(topsnps_PAGER$topsnp, split=":"), `[`, 2)

# Initialize an empty dataframe to store results
out <- data.frame()

# Loop over each SNP to calculate LD intervals
for(i in 1:nrow(topsnps_PAGER)) {
  trait <- topsnps_PAGER$trait[i]
  topsnp <- topsnps_PAGER$topsnp[i]
  chr <- topsnps_PAGER$chr[i]
  pos <- topsnps_PAGER$pos[i]
  
  # Run PLINK to compute LD, outputting results to a temporary file
  system(paste0("./plink --bfile ",exp_dir,"/P50_round2_LD_pruned_3473 --chr ",gsub("chr", "", chr)," --nonfounders --r2  --ld-snp ",topsnp," --ld-window 100000 --ld-window-kb 6000 --ld-window-r2 0 --out ",exp_dir,"/LZ_temp"),ignore.stdout = T,ignore.stderr = T,wait = T)
  
  # Calculate LD interval start and stop positions
  test <- LD_start_stop(paste0(exp_dir, "/LZ_temp.ld"), trait = trait, chr = chr, pos = pos)
  out <- rbind(out, test) # Append results to the output dataframe
}

# Bring LD interval information into the dataframe
topsnps_PAGER$LD_start <- out$LD_interval_start
topsnps_PAGER$LD_stop <- out$LD_interval_stop

# Write updated file to disk
write.csv(topsnps_PAGER, "PAGER_QTLs_LD.csv", row.names = FALSE)
