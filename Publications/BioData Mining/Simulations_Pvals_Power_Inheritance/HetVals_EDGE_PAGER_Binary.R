### R
### PJF
### 3/14/24

### This script will calculate average heterozygote values between EDGE and PAGER for comparison purposes. Discrete Phenotype

## Discrete
setwd("/path/to/folder")

# Load necessary packages
library(dplyr)
library(tidyr)
library(broom)
library(tidyverse)

# Function to extract the first string before the first underscore from the filename
extract_name <- function(filename) {
  parts <- unlist(strsplit(filename, "_", fixed = TRUE))
  # Return the first segment before any underscore or the filename itself if no underscore exists
  return(parts[1])
}

# Get a list of all CSV files in the current working directory
csv_files <- list.files(pattern="*.csv")

# Loop through each file
for(file in csv_files) {
  # Extract a valid R object name from the filename
  dataframe_name <- extract_name(file)
  
  # Read the CSV file into a temporary variable
  temp_dataframe <- read_csv(file)
  
  # Assign the dataframe to a dynamically named variable in the global environment
  assign(dataframe_name, temp_dataframe, envir = .GlobalEnv)
  
  # Remove the temporary dataframe variable to avoid clutter
  rm(temp_dataframe)
}

# Add in Inherent p-value for additive so we can Rbind
additive$INHERENT_pvalue <- NA
additive <- additive[, c(1,2,3,4,5,6,7,8,9,10,15,11,12,13,14)]

# Combine dataframes
dataframes <- list(
  additive = additive,
  dominant = dominant,
  heterosis = heterosis,
  overdominant = overdominant,
  recessive = recessive,
  subadditive = subadditive,
  superadditive = superadditive,
  underdominant = underdominant
)

# Combine the dataframes with an additional column indicating the original dataframe name
combined_discrete <- bind_rows(dataframes, .id = "Inheritance")

## All Data
# Calculate averages for each inheritance Model
averages <- combined_discrete %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_EDGE = mean(EDGE_VALUE, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_VALUE_1, na.rm = TRUE),
    Median_of_EDGE = median(EDGE_VALUE, na.rm = TRUE),
    Median_of_PAGER = median(PAGER_VALUE_1, na.rm = TRUE),
    Var_of_EDGE = var(EDGE_VALUE, na.rm = TRUE),
    Var_of_PAGER = var(PAGER_VALUE_1, na.rm = TRUE),
    SD_of_EDGE = sd(EDGE_VALUE, na.rm = TRUE),
    SD_of_PAGER = sd(PAGER_VALUE_1, na.rm = TRUE),
    Min_of_EDGE = min(EDGE_VALUE, na.rm = TRUE),
    Min_of_PAGER = min(PAGER_VALUE_1, na.rm = TRUE),
    Max_of_EDGE = max(EDGE_VALUE, na.rm = TRUE),
    Max_of_PAGER = max(PAGER_VALUE_1, na.rm = TRUE),
    Range_of_EDGE = max(EDGE_VALUE, na.rm = TRUE) - min(EDGE_VALUE, na.rm = TRUE),
    Range_of_PAGER = max(PAGER_VALUE_1, na.rm = TRUE) - min(PAGER_VALUE_1, na.rm = TRUE),
    IQR_of_EDGE = IQR(EDGE_VALUE, na.rm = TRUE),
    IQR_of_PAGER = IQR(PAGER_VALUE_1, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data <- pivot_longer(combined_discrete, cols = c("EDGE_VALUE", "PAGER_VALUE_1"),
                          names_to = "Encoding", values_to = "Het_Value")

# Perform Wilcoxon Rank SUM Tests grouped by Inheritance_Model and PValue_Type
Wilcox_test_results <- long_data %>%
  group_by(Inheritance) %>%
  do(tidy(wilcox.test(Het_Value ~ Encoding, data = .)))

# Write to disk
setwd("/path/to/folder")
write.csv(averages, "Binary_AvgHet_All.csv", row.names = FALSE)
write.csv(Wilcox_test_results, "Binary_WilcoxHet_All.csv", row.names = FALSE)

## Within Inheritance Model for MAF = 0.1
Discrete_filtered1 <- combined_discrete %>% filter(MAFA == 0.1)

# Calculate averages for each inheritance Model
averages1 <- Discrete_filtered1 %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_EDGE = mean(EDGE_VALUE, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_VALUE_1, na.rm = TRUE),
    Median_of_EDGE = median(EDGE_VALUE, na.rm = TRUE),
    Median_of_PAGER = median(PAGER_VALUE_1, na.rm = TRUE),
    Var_of_EDGE = var(EDGE_VALUE, na.rm = TRUE),
    Var_of_PAGER = var(PAGER_VALUE_1, na.rm = TRUE),
    SD_of_EDGE = sd(EDGE_VALUE, na.rm = TRUE),
    SD_of_PAGER = sd(PAGER_VALUE_1, na.rm = TRUE),
    Min_of_EDGE = min(EDGE_VALUE, na.rm = TRUE),
    Min_of_PAGER = min(PAGER_VALUE_1, na.rm = TRUE),
    Max_of_EDGE = max(EDGE_VALUE, na.rm = TRUE),
    Max_of_PAGER = max(PAGER_VALUE_1, na.rm = TRUE),
    Range_of_EDGE = max(EDGE_VALUE, na.rm = TRUE) - min(EDGE_VALUE, na.rm = TRUE),
    Range_of_PAGER = max(PAGER_VALUE_1, na.rm = TRUE) - min(PAGER_VALUE_1, na.rm = TRUE),
    IQR_of_EDGE = IQR(EDGE_VALUE, na.rm = TRUE),
    IQR_of_PAGER = IQR(PAGER_VALUE_1, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data1 <- pivot_longer(Discrete_filtered1, cols = c("EDGE_VALUE", "PAGER_VALUE_1"),
                          names_to = "Encoding", values_to = "Het_Value")

# Perform Wilcoxon Rank SUM Tests grouped by Inheritance_Model and PValue_Type
Wilcox_test_results1 <- long_data1 %>%
  group_by(Inheritance) %>%
  do(tidy(wilcox.test(Het_Value ~ Encoding, data = .)))

# Write to disk
setwd("/path/to/folder")
write.csv(averages1, "Binary_AvgHet_MAF01.csv", row.names = FALSE)
write.csv(Wilcox_test_results1, "Binary_WilcoxHet_MAF01.csv", row.names = FALSE)

## Within Inheritance Model for MAF = 0.1 and PEN_DIFF = 0.1
Discrete_filtered2 <- Discrete_filtered1  %>% filter(PEN_DIFF == 0.1)

# Calculate averages for each inheritance Model
averages2 <- Discrete_filtered2 %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_EDGE = mean(EDGE_VALUE, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_VALUE_1, na.rm = TRUE),
    Median_of_EDGE = median(EDGE_VALUE, na.rm = TRUE),
    Median_of_PAGER = median(PAGER_VALUE_1, na.rm = TRUE),
    Var_of_EDGE = var(EDGE_VALUE, na.rm = TRUE),
    Var_of_PAGER = var(PAGER_VALUE_1, na.rm = TRUE),
    SD_of_EDGE = sd(EDGE_VALUE, na.rm = TRUE),
    SD_of_PAGER = sd(PAGER_VALUE_1, na.rm = TRUE),
    Min_of_EDGE = min(EDGE_VALUE, na.rm = TRUE),
    Min_of_PAGER = min(PAGER_VALUE_1, na.rm = TRUE),
    Max_of_EDGE = max(EDGE_VALUE, na.rm = TRUE),
    Max_of_PAGER = max(PAGER_VALUE_1, na.rm = TRUE),
    Range_of_EDGE = max(EDGE_VALUE, na.rm = TRUE) - min(EDGE_VALUE, na.rm = TRUE),
    Range_of_PAGER = max(PAGER_VALUE_1, na.rm = TRUE) - min(PAGER_VALUE_1, na.rm = TRUE),
    IQR_of_EDGE = IQR(EDGE_VALUE, na.rm = TRUE),
    IQR_of_PAGER = IQR(PAGER_VALUE_1, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data2 <- pivot_longer(Discrete_filtered2, cols = c("EDGE_VALUE", "PAGER_VALUE_1"),
                           names_to = "Encoding", values_to = "Het_Value")

# Perform Wilcoxon Rank SUM Tests grouped by Inheritance_Model and PValue_Type
Wilcox_test_results2 <- long_data2 %>%
  group_by(Inheritance) %>%
  do(tidy(wilcox.test(Het_Value ~ Encoding, data = .)))

# Write to disk
setwd("/path/to/folder")
write.csv(averages2, "Binary_AvgHet_MAF01_PD01.csv", row.names = FALSE)
write.csv(Wilcox_test_results2, "Binary_WilcoxHet_MAF01_PD01.csv", row.names = FALSE)
