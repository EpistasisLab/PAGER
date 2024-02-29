setwd("/Users/fredap/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/FPR/null simulation")

library(tidyverse)

# Step 1: Read all CSV files and combine them
file_names <- list.files(path = ".", pattern = "null_effect_simulation_results.*\\.csv")
combined_df <- file_names %>% 
  lapply(read.csv) %>% 
  bind_rows()

write.csv(combined_df, "combined_df.csv", row.names = FALSE)

tests_per_group <- 8000 # Assuming we aggregate across C/C

# Step 2: Filter and Aggregate Data
result_df <- combined_df %>%
  group_by(MAFA, Num_Samples) %>%
  summarize(
    Num_additive_sig = sum(additive_pvalue < 0.05, na.rm = TRUE),
    Num_EDGE_sig = sum(EDGE_pvalue < 0.05, na.rm = TRUE),
    Num_PAGER_sig = sum(PAGER_pvalue < 0.05, na.rm = TRUE)
  ) %>%
  mutate(
    FPR_additive = Num_additive_sig / tests_per_group,
    FPR_EDGE = Num_EDGE_sig / tests_per_group,
    FPR_PAGER = Num_PAGER_sig / tests_per_group
  )

write.csv(result_df, "FPR_per_Encoding.csv", row.names = FALSE)