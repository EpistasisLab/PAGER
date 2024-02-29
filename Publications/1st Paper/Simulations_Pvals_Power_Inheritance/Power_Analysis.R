# Power Analysis 

library(dplyr)

# Additive encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_additive_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, Case_Control_Ratio, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "additive_power_results.csv", row.names = FALSE)

#################################################

# Dominant encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_dominant_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, Case_Control_Ratio, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_pwer = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "dominant_power_results.csv", row.names = FALSE)


#################################################

# Recessive encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_recessive_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, Case_Control_Ratio, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_pwer = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "recessive_power_results.csv", row.names = FALSE)


#################################################

# Heterosis encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_heterosis_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, Case_Control_Ratio, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_pwer = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "heterosis_power_results.csv", row.names = FALSE)


#################################################

# Superadditive encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_superadditive_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, Case_Control_Ratio, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_pwer = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "superadditive_power_results.csv", row.names = FALSE)


#################################################

# Subadditive encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_subadditive_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, Case_Control_Ratio, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_pwer = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "subadditive_power_results.csv", row.names = FALSE)


#################################################

# Overdominant encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_overdominant_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, Case_Control_Ratio, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_pwer = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "overdominant_power_results.csv", row.names = FALSE)


#################################################

# Underdominant encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/clarite_underdominant_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, Case_Control_Ratio, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_power = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "underdominant_power_results.csv", row.names = FALSE)


#############################################################################################################################

############# POWER ANALYSIS NON-DISCRETE/CONTINUOUS ###############

# Additive encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_additive_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/additive_power_results.csv", row.names = FALSE)

#################################################

# Dominant encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_dominant_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_power = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/dominant_power_results.csv", row.names = FALSE)


#################################################

# Recessive encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_recessive_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_power = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/recessive_power_results.csv", row.names = FALSE)


#################################################

# Heterosis encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_heterosis_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_power = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/heterosis_power_results.csv", row.names = FALSE)


#################################################

# Superadditive encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_superadditive_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_power = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/superadditive_power_results.csv", row.names = FALSE)


#################################################

# Subadditive encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_subadditive_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_power = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/subadditive_power_results.csv", row.names = FALSE)


#################################################

# Overdominant encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_overdominant_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_power = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/overdominant_power_results.csv", row.names = FALSE)


#################################################

# Underdominant encoding
# read the file
data <- read.csv("/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/clarite_cont_underdominant_simulation_results.csv")

# Assuming 'data' is your dataframe
result <- data %>%
  group_by(Num_Samples, MAFA, PEN_DIFF) %>%
  summarize(
    Original_Less_T = sum(Original_pvalue < 5e-8),
    Original_Greater_F = sum(Original_pvalue >= 5e-8),
    INHERENT_Less_T = sum(INHERENT_pvalue < 5e-8),
    INHERENT_Greater_F = sum(INHERENT_pvalue >= 5e-8),
    EDGE_Less_T = sum(EDGE_pvalue < 5e-8),
    EDGE_Greater_F = sum(EDGE_pvalue >= 5e-8),
    PAGER_Less_T = sum(PAGER_pvalue < 5e-8),
    PAGER_Greater_F = sum(PAGER_pvalue >= 5e-8)
  ) %>%
  ungroup()

# Assuming 'result' is your dataframe with the counts

result <- result %>%
  mutate(
    Original_power = Original_Less_T / 1000,
    INHERENT_power = INHERENT_Less_T /1000,
    EDGE_power = EDGE_Less_T / 1000,
    PAGER_power = PAGER_Less_T / 1000
  )

# Save the results to a CSV file
write.csv(result, "/Users/ghosha/Library/CloudStorage/Box-Box/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_non_discrete/underdominant_power_results.csv", row.names = FALSE)


###############################################################################################
