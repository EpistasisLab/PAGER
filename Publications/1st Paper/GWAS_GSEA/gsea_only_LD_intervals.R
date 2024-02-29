# This is the script to conduct GSEA on the two LD intervals of interest

library(RJSONIO)
library(tidyr)
library(dplyr)
library(stringr)
library(cowplot)
library(data.table)
library(BiocManager)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(ggplot2)
library(org.Rn.eg.db)

#set directories
mainDir <- "/common/compbiomed-dsn/projects/Philip_Freda_GWAS_20231024/analysis/gsea/"
mainDir <- "/Users/bhandaryp/Documents/philip_freda_gwas/"
setwd(mainDir)

# experiment 1 - LD interval 1
# chr18.26640423 - chr18.27355039
# define the chromosome, and upstream and downstream basepairs from LD interval
chr <- 18
upstream <- 26640423
downstream <- 27355039

url <- paste("https://rest.rgd.mcw.edu/rgdws/genes/", chr, "/", upstream, "/" , downstream, "/360", sep = "") # Use the URL from the RGD database to obtain genes in that interval
print(url)
holder_exp1 <- as.data.frame(lapply(url, jsonlite::fromJSON, flatten= TRUE)) # Use lapply to get the corresponding genes into a dataframe
holder_exp1$locus_name <- "chr18.26640423-chr18.27355039" # Also add the LD intervals and QTL information if you have any
holder_exp1$QTL <- "FALSE"
gene_list <- holder_exp1$symbol # Add the gene symbol from the dataframe to a list that can be used

hs <- org.Rn.eg.db # Get the Rat database annotations
my.symbols <- gene_list # Save your gene symbols as my.symbols 
Entrez <- select(hs,
                 keys = my.symbols,
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = "SYMBOL") # Use select from AnnotationDbi to get the Entrez gene IDs for the genes

Entrez <- na.omit(Entrez) # remove any rows that contain NA's
write.table(Entrez, "LD_intervals_experiment1_Enrichr.csv", row.names = FALSE) # write out Entrez information for Enrichr upload
gene_list <- Entrez$ENTREZID # save this column into gene_list

# KEGG Pathway Analysis
onlyKEGG <- enrichKEGG(gene = gene_list, organism = 'rno') # Use enrichKEGG function to get the KEGG pathway enrichment
onlyKEGG@result <- onlyKEGG@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the onlyKEGG object to get only those pathways where the pvalue and qvalue is less than 0.05
write.csv(onlyKEGG@result, "LD_intervals_experiment1_KEGG_only.csv", row.names = FALSE) # Write to disk

ego_CC <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont CC to peform GO enrichment for CC ontology

ego_CC@result <- ego_CC@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_CC object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_CC@result, "LD_intervals_experiment1_CC_only.csv", row.names = FALSE) # Write to disk 

ego_BP <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont BP to peform GO enrichment for BP ontology

ego_BP@result <- ego_BP@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_BP object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_BP@result, "LD_intervals_experiment1_BP_only.csv", row.names = FALSE) # Write to disk

ego_MF <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont MF to peform GO enrichment for BP ontology

ego_MF@result <- ego_MF@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_MF object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_MF@result, "LD_intervals_experiment1_MF_only.csv", row.names = FALSE) # Write to disk

# experiment 2 - LD interval 2
# chr8.81203483 - chr8.81514453
chr <- 8
upstream <- 81203483
downstream <- 81514453


url <- paste("https://rest.rgd.mcw.edu/rgdws/genes/", chr, "/", upstream, "/" , downstream, "/360", sep = "")
print(url)
holder_exp2 <- as.data.frame(lapply(url, jsonlite::fromJSON, flatten= TRUE))
holder_exp2$locus_name <- "chr8.81203483 - chr8.81514453"
holder_exp2$QTL <- "FALSE"
gene_list <- holder_exp2$symbol


hs <- org.Rn.eg.db # Get the Rat database annotations
my.symbols <- gene_list # Save your gene symbols as my.symbols 
Entrez <- select(hs,
                 keys = my.symbols,
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = "SYMBOL") # Use select from AnnotationDbi to get the Entrez gene IDs for the genes

Entrez <- na.omit(Entrez) # remove any rows that contain NA's
write.table(Entrez, "LD_intervals_experiment2_Enrichr.csv", row.names = FALSE) # write out Entrez information for Enrichr upload
gene_list <- Entrez$ENTREZID # save this column into gene_list

# KEGG Pathway Analysis
onlyKEGG <- enrichKEGG(gene = gene_list, organism = 'rno') # Use enrichKEGG function to get the KEGG pathway enrichment
onlyKEGG@result <- onlyKEGG@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the onlyKEGG object to get only those pathways where the pvalue and qvalue is less than 0.05
write.csv(onlyKEGG@result, "LD_intervals_experiment2_KEGG_only.csv", row.names = FALSE) # Write to disk

ego_CC <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont CC to peform GO enrichment for CC ontology

ego_CC@result <- ego_CC@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_CC object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_CC@result, "LD_intervals_experiment2_CC_only.csv", row.names = FALSE) # Write to disk 

ego_BP <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont BP to peform GO enrichment for BP ontology

ego_BP@result <- ego_BP@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_BP object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_BP@result, "LD_intervals_experiment2_BP_only.csv", row.names = FALSE) # Write to disk

ego_MF <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont MF to peform GO enrichment for BP ontology

ego_MF@result <- ego_MF@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_MF object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_MF@result, "LD_intervals_experiment2_MF_only.csv", row.names = FALSE) # Write to disk


# experiment 3 - LD interval 1 + LD interval 2
holder_exp3 <- rbind(holder_exp1, holder_exp2)
gene_list <- holder_exp3$symbol

hs <- org.Rn.eg.db # Get the Rat database annotations
my.symbols <- gene_list # Save your gene symbols as my.symbols 
Entrez <- select(hs,
                 keys = my.symbols,
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = "SYMBOL") # Use select from AnnotationDbi to get the Entrez gene IDs for the genes

Entrez <- na.omit(Entrez) # remove any rows that contain NA's
write.table(Entrez, "LD_intervals_experiment3_Enrichr.csv", row.names = FALSE) # write out Entrez information for Enrichr upload
gene_list <- Entrez$ENTREZID # save this column into gene_list

# KEGG Pathway Analysis
onlyKEGG <- enrichKEGG(gene = gene_list, organism = 'rno') # Use enrichKEGG function to get the KEGG pathway enrichment
onlyKEGG@result <- onlyKEGG@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the onlyKEGG object to get only those pathways where the pvalue and qvalue is less than 0.05
write.csv(onlyKEGG@result, "LD_intervals_experiment3_KEGG_only.csv", row.names = FALSE) # Write to disk

ego_CC <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont CC to peform GO enrichment for CC ontology

ego_CC@result <- ego_CC@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_CC object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_CC@result, "LD_intervals_experiment3_CC_only.csv", row.names = FALSE) # Write to disk 

ego_BP <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont BP to peform GO enrichment for BP ontology

ego_BP@result <- ego_BP@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_BP object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_BP@result, "LD_intervals_experiment3_BP_only.csv", row.names = FALSE) # Write to disk

ego_MF <- enrichGO(gene          = gene_list,
                   OrgDb         = org.Rn.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "fdr",
                   readable      = TRUE) # Use enrichGO function and use ont MF to peform GO enrichment for BP ontology

ego_MF@result <- ego_MF@result %>% filter(pvalue < 0.05, qvalue < 0.05) # Filter the result dataframe inside of the ego_MF object to get only those GO terms where the pvalue and qvalue is less than 0.05
write.csv(ego_MF@result, "LD_intervals_experiment3_MF_only.csv", row.names = FALSE) # Write to disk