# Power Analysis 

library(dplyr)

################## FOR DISCRETE PHENOTYPE FILES ################

# Additive encoding
# read the file
data <- read.csv("/path/to/read/additive_simulation_results.csv")

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
write.csv(result, "/path/to/save/additive_power_results.csv", row.names = FALSE)

#################################################

# Dominant encoding
# read the file
data <- read.csv("/path/to/read/dominant_simulation_results.csv")

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
write.csv(result, "/path/to/save/dominant_power_results.csv", row.names = FALSE)


#################################################

# Recessive encoding
# read the file
data <- read.csv("/path/to/read/recessive_simulation_results.csv")

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
write.csv(result, "/path/to/save/recessive_power_results.csv", row.names = FALSE)


#################################################

# Heterosis encoding
# read the file
data <- read.csv("/path/to/read/heterosis_simulation_results.csv")

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
write.csv(result, "/path/to/save/heterosis_power_results.csv", row.names = FALSE)


#################################################

# Superadditive encoding
# read the file
data <- read.csv("/path/to/read/superadditive_simulation_results.csv")

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
write.csv(result, "/path/to/save/superadditive_power_results.csv", row.names = FALSE)


#################################################

# Subadditive encoding
# read the file
data <- read.csv("/path/to/read/subadditive_simulation_results.csv")

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
write.csv(result, "/path/to/save/subadditive_power_results.csv", row.names = FALSE)


#################################################

# Overdominant encoding
# read the file
data <- read.csv("/path/to/read/overdominant_simulation_results.csv")

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
write.csv(result, "/path/to/save/overdominant_power_results.csv", row.names = FALSE)


#################################################

# Underdominant encoding
# read the file
data <- read.csv("/path/to/read/underdominant_simulation_results.csv")

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
write.csv(result, "/path/to/save/underdominant_power_results.csv", row.names = FALSE)


#############################################################################################################################

############# FOR NON-DISCRETE/CONTINUOUS FILES ###############

# Additive encoding
# read the file
data <- read.csv("/path/to/read/cont_additive_simulation_results.csv")

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
write.csv(result, "/path/to/save/additive_power_results.csv", row.names = FALSE)

#################################################

# Dominant encoding
# read the file
data <- read.csv("/path/to/read/cont_dominant_simulation_results.csv")

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
write.csv(result, "/path/to/save/dominant_power_results.csv", row.names = FALSE)


#################################################

# Recessive encoding
# read the file
data <- read.csv("/path/to/read/cont_recessive_simulation_results.csv")

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
write.csv(result, "/path/to/save/recessive_power_results.csv", row.names = FALSE)


#################################################

# Heterosis encoding
# read the file
data <- read.csv("/path/to/read/cont_heterosis_simulation_results.csv")

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
write.csv(result, "/path/to/save/heterosis_power_results.csv", row.names = FALSE)


#################################################

# Superadditive encoding
# read the file
data <- read.csv("/path/to/read/cont_superadditive_simulation_results.csv")

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
write.csv(result, "/path/to/save/superadditive_power_results.csv", row.names = FALSE)


#################################################

# Subadditive encoding
# read the file
data <- read.csv("/path/to/read/cont_subadditive_simulation_results.csv")

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
write.csv(result, "/path/to/save/subadditive_power_results.csv", row.names = FALSE)


#################################################

# Overdominant encoding
# read the file
data <- read.csv("/path/to/read/cont_overdominant_simulation_results.csv")

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
write.csv(result, "/path/to/save/overdominant_power_results.csv", row.names = FALSE)


#################################################

# Underdominant encoding
# read the file
data <- read.csv("/path/to/read/cont_underdominant_simulation_results.csv")

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
write.csv(result, "/path/to/save/underdominant_power_results.csv", row.names = FALSE)


###############################################################################################
