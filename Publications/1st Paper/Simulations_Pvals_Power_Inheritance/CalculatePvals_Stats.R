library(ggplot2)
library(dplyr)

## Continuous

# Import files
add <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_additive_simulation_results.csv")
dom <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_dominant_simulation_results.csv")
rec <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_recessive_simulation_results.csv")
subadd <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_subadditive_simulation_results.csv")
supadd <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_superadditive_simulation_results.csv")
het <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_heterosis_simulation_results.csv")
und <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_underdominant_simulation_results.csv")
ovd <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_overdominant_simulation_results.csv")

calculate_stats <- function(column) {
  n <- sum(!is.na(column)) # Count of non-NA values for sample size
  stats <- c(
    mean = mean(column, na.rm = TRUE),
    median = median(column, na.rm = TRUE),
    sd = sd(column, na.rm = TRUE), # Standard deviation
    var = var(column, na.rm = TRUE), # Variance
    sem = sd(column, na.rm = TRUE) / sqrt(n) # Standard error of the mean
  )
  return(stats)
}

# Additive IM
Original_pvalue_stats_add <- calculate_stats(add$Original_pvalue)
EDGE_pvalue_stats_add <- calculate_stats(add$EDGE_pvalue)
PAGER_pvalue_stats_add <- calculate_stats(add$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_add <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_add),
  Inherent_pvalue = -log10(Original_pvalue_stats_add),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_add),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_add),
  row.names = NULL # This ensures no row names are set
)

stats_df_add$InheritanceModel <- 'Additive'

# Dominant IM
Original_pvalue_stats_dom <- calculate_stats(dom$Original_pvalue)
Inherent_pvalue_stats_dom <- calculate_stats(dom$INHERENT_pvalue)
EDGE_pvalue_stats_dom <- calculate_stats(dom$EDGE_pvalue)
PAGER_pvalue_stats_dom <- calculate_stats(dom$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_dom <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_dom),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_dom),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_dom),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_dom),
  row.names = NULL # This ensures no row names are set
)

stats_df_dom$InheritanceModel <- 'Dominant'

# Recessive IM
Original_pvalue_stats_rec <- calculate_stats(rec$Original_pvalue)
Inherent_pvalue_stats_rec <- calculate_stats(rec$INHERENT_pvalue)
EDGE_pvalue_stats_rec <- calculate_stats(rec$EDGE_pvalue)
PAGER_pvalue_stats_rec <- calculate_stats(rec$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_rec <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_rec),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_rec),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_rec),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_rec),
  row.names = NULL # This ensures no row names are set
)

stats_df_rec$InheritanceModel <- 'Recessive'

# Subadditive IM
Original_pvalue_stats_subadd <- calculate_stats(subadd$Original_pvalue)
Inherent_pvalue_stats_subadd <- calculate_stats(subadd$INHERENT_pvalue)
EDGE_pvalue_stats_subadd <- calculate_stats(subadd$EDGE_pvalue)
PAGER_pvalue_stats_subadd <- calculate_stats(subadd$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_subadd <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_subadd),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_subadd),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_subadd),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_subadd),
  row.names = NULL # This ensures no row names are set
)

stats_df_subadd$InheritanceModel <- 'Subadditive'

# Superadditive IM
Original_pvalue_stats_supadd <- calculate_stats(supadd$Original_pvalue)
Inherent_pvalue_stats_supadd <- calculate_stats(supadd$INHERENT_pvalue)
EDGE_pvalue_stats_supadd <- calculate_stats(supadd$EDGE_pvalue)
PAGER_pvalue_stats_supadd <- calculate_stats(supadd$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_supadd <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_supadd),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_supadd),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_supadd),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_supadd),
  row.names = NULL # This ensures no row names are set
)

stats_df_supadd$InheritanceModel <- 'Superadditive'

# Heterosis IM
Original_pvalue_stats_het <- calculate_stats(het$Original_pvalue)
Inherent_pvalue_stats_het <- calculate_stats(het$INHERENT_pvalue)
EDGE_pvalue_stats_het <- calculate_stats(het$EDGE_pvalue)
PAGER_pvalue_stats_het <- calculate_stats(het$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_het <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_het),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_het),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_het),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_het),
  row.names = NULL # This ensures no row names are set
)

stats_df_het$InheritanceModel <- 'Heterosis'

# Underdominant IM
Original_pvalue_stats_und <- calculate_stats(und$Original_pvalue)
Inherent_pvalue_stats_und <- calculate_stats(und$INHERENT_pvalue)
EDGE_pvalue_stats_und <- calculate_stats(und$EDGE_pvalue)
PAGER_pvalue_stats_und <- calculate_stats(und$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_und <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_und),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_und),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_und),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_und),
  row.names = NULL # This ensures no row names are set
)

stats_df_und$InheritanceModel <- 'Underdominant'

# Overdominant IM
Original_pvalue_stats_ovd <- calculate_stats(ovd$Original_pvalue)
Inherent_pvalue_stats_ovd <- calculate_stats(ovd$INHERENT_pvalue)
EDGE_pvalue_stats_ovd <- calculate_stats(ovd$EDGE_pvalue)
PAGER_pvalue_stats_ovd <- calculate_stats(ovd$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_ovd <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_ovd),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_ovd),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_ovd),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_ovd),
  row.names = NULL # This ensures no row names are set
)

stats_df_ovd$InheritanceModel <- 'Overdominant'

# rbind all together
combined_stats_df <- rbind(stats_df_add, stats_df_subadd, stats_df_supadd, stats_df_dom, stats_df_rec, stats_df_het, stats_df_und, stats_df_ovd)

# write to disk
write.csv(combined_stats_df, "cont_pvals.csv", row.names = FALSE)

## Discrete

# Import files
add <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_additive_simulation_results.csv")
dom <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_dominant_simulation_results.csv")
rec <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_recessive_simulation_results.csv")
subadd <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_subadditive_simulation_results.csv")
supadd <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_superadditive_simulation_results.csv")
het <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_heterosis_simulation_results.csv")
und <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_underdominant_simulation_results.csv")
ovd <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_overdominant_simulation_results.csv")

calculate_stats <- function(column) {
  n <- sum(!is.na(column)) # Count of non-NA values for sample size
  stats <- c(
    mean = mean(column, na.rm = TRUE),
    median = median(column, na.rm = TRUE),
    sd = sd(column, na.rm = TRUE), # Standard deviation
    var = var(column, na.rm = TRUE), # Variance
    sem = sd(column, na.rm = TRUE) / sqrt(n) # Standard error of the mean
  )
  return(stats)
}

# Additive IM
Original_pvalue_stats_add <- calculate_stats(add$Original_pvalue)
EDGE_pvalue_stats_add <- calculate_stats(add$EDGE_pvalue)
PAGER_pvalue_stats_add <- calculate_stats(add$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_add <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_add),
  Inherent_pvalue = -log10(Original_pvalue_stats_add),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_add),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_add),
  row.names = NULL # This ensures no row names are set
)

stats_df_add$InheritanceModel <- 'Additive'

# Dominant IM
Original_pvalue_stats_dom <- calculate_stats(dom$Original_pvalue)
Inherent_pvalue_stats_dom <- calculate_stats(dom$INHERENT_pvalue)
EDGE_pvalue_stats_dom <- calculate_stats(dom$EDGE_pvalue)
PAGER_pvalue_stats_dom <- calculate_stats(dom$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_dom <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_dom),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_dom),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_dom),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_dom),
  row.names = NULL # This ensures no row names are set
)

stats_df_dom$InheritanceModel <- 'Dominant'

# Recessive IM
Original_pvalue_stats_rec <- calculate_stats(rec$Original_pvalue)
Inherent_pvalue_stats_rec <- calculate_stats(rec$INHERENT_pvalue)
EDGE_pvalue_stats_rec <- calculate_stats(rec$EDGE_pvalue)
PAGER_pvalue_stats_rec <- calculate_stats(rec$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_rec <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_rec),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_rec),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_rec),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_rec),
  row.names = NULL # This ensures no row names are set
)

stats_df_rec$InheritanceModel <- 'Recessive'

# Subadditive IM
Original_pvalue_stats_subadd <- calculate_stats(subadd$Original_pvalue)
Inherent_pvalue_stats_subadd <- calculate_stats(subadd$INHERENT_pvalue)
EDGE_pvalue_stats_subadd <- calculate_stats(subadd$EDGE_pvalue)
PAGER_pvalue_stats_subadd <- calculate_stats(subadd$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_subadd <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_subadd),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_subadd),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_subadd),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_subadd),
  row.names = NULL # This ensures no row names are set
)

stats_df_subadd$InheritanceModel <- 'Subadditive'

# Superadditive IM
Original_pvalue_stats_supadd <- calculate_stats(supadd$Original_pvalue)
Inherent_pvalue_stats_supadd <- calculate_stats(supadd$INHERENT_pvalue)
EDGE_pvalue_stats_supadd <- calculate_stats(supadd$EDGE_pvalue)
PAGER_pvalue_stats_supadd <- calculate_stats(supadd$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_supadd <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_supadd),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_supadd),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_supadd),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_supadd),
  row.names = NULL # This ensures no row names are set
)

stats_df_supadd$InheritanceModel <- 'Superadditive'

# Heterosis IM
Original_pvalue_stats_het <- calculate_stats(het$Original_pvalue)
Inherent_pvalue_stats_het <- calculate_stats(het$INHERENT_pvalue)
EDGE_pvalue_stats_het <- calculate_stats(het$EDGE_pvalue)
PAGER_pvalue_stats_het <- calculate_stats(het$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_het <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_het),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_het),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_het),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_het),
  row.names = NULL # This ensures no row names are set
)

stats_df_het$InheritanceModel <- 'Heterosis'

# Underdominant IM
Original_pvalue_stats_und <- calculate_stats(und$Original_pvalue)
Inherent_pvalue_stats_und <- calculate_stats(und$INHERENT_pvalue)
EDGE_pvalue_stats_und <- calculate_stats(und$EDGE_pvalue)
PAGER_pvalue_stats_und <- calculate_stats(und$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_und <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_und),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_und),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_und),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_und),
  row.names = NULL # This ensures no row names are set
)

stats_df_und$InheritanceModel <- 'Underdominant'

# Overdominant IM
Original_pvalue_stats_ovd <- calculate_stats(ovd$Original_pvalue)
Inherent_pvalue_stats_ovd <- calculate_stats(ovd$INHERENT_pvalue)
EDGE_pvalue_stats_ovd <- calculate_stats(ovd$EDGE_pvalue)
PAGER_pvalue_stats_ovd <- calculate_stats(ovd$PAGER_pvalue)

# Combine the stats into a new dataframe
stats_df_ovd <- data.frame(
  Statistic = c('mean', 'median', 'sd', 'var', 'sem'),
  Original_pvalue = -log10(Original_pvalue_stats_ovd),
  Inherent_pvalue = -log10(Inherent_pvalue_stats_ovd),
  EDGE_pvalue = -log10(EDGE_pvalue_stats_ovd),
  PAGER_pvalue = -log10(PAGER_pvalue_stats_ovd),
  row.names = NULL # This ensures no row names are set
)

stats_df_ovd$InheritanceModel <- 'Overdominant'

# rbind all together
combined_stats_df <- rbind(stats_df_add, stats_df_subadd, stats_df_supadd, stats_df_dom, stats_df_rec, stats_df_het, stats_df_und, stats_df_ovd)

# write to disk
write.csv(combined_stats_df, "disc_pvals.csv", row.names = FALSE)
