### R
### PJF
### 3/14/24

### This script will show average p-values for a host of certain combinations within inheritance models

## Discrete
setwd("/Users/fredap/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete")

# Load necessary packages
library(dplyr)
library(tidyr)
library(broom)
library(tidyverse)

# Function to extract the desired pattern from the filename
extract_name <- function(filename) {
  parts <- unlist(strsplit(filename, "_"))
  if (length(parts) > 1) {
    # Return the part of the filename between the first set of underscores
    return(parts[2])
  } else {
    # If there's no underscore, return the filename without the extension
    return(tools::file_path_sans_ext(filename))
  }
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

# Multiply EDGE pvalue by 2 to control for multiple testing
combined_discrete$EDGE_pvalue <- combined_discrete$EDGE_pvalue*2

## All Data
# Calculate average p-vals for each inheritance Model
averages_filtered <- combined_discrete %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_Original = mean(Original_pvalue, na.rm = TRUE),
    Avg_of_INHERENT = mean(INHERENT_pvalue, na.rm = TRUE),
    Avg_of_EDGE = mean(EDGE_pvalue, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_pvalue, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data <- pivot_longer(combined_discrete, cols = c("Original_pvalue", "INHERENT_pvalue", "EDGE_pvalue", "PAGER_pvalue"),
                           names_to = "PValue_Type", values_to = "PValue") %>%
  filter(!(Inheritance == "additive" & PValue_Type == "INHERENT_pvalue"))

# Perform ANOVA grouped by Inheritance_Model and PValue_Type
anova_results <- long_data %>%
  group_by(Inheritance) %>%
  do(tidy(aov(PValue ~ PValue_Type, data = .)))

# Perform Tukey's HSD test based on the corrected ANOVA setup
tukey_results <- long_data %>%
  group_by(Inheritance) %>%
  do({
    model <- aov(PValue ~ PValue_Type, data = .)
    tukey <- TukeyHSD(model)
    tidy(tukey)
  })

# Write to disk
write.csv(averages_filtered, "All_Binary_AvgP.csv", row.names = FALSE)
write.csv(anova_results, "All_Binary_ANOVA.csv", row.names = FALSE)
write.csv(tukey_results, "All_Binary_Tukey.csv", row.names = FALSE)

## Within Inheritance Model for MAF = 0.1
Discrete_filtered1 <- combined_discrete %>% filter(MAFA == 0.1)

# Calculate average p-vals for each inheritance Model
averages_filtered1 <- Discrete_filtered1 %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_Original = mean(Original_pvalue, na.rm = TRUE),
    Avg_of_INHERENT = mean(INHERENT_pvalue, na.rm = TRUE),
    Avg_of_EDGE = mean(EDGE_pvalue, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_pvalue, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data1 <- pivot_longer(Discrete_filtered1, cols = c("Original_pvalue", "INHERENT_pvalue", "EDGE_pvalue", "PAGER_pvalue"),
                          names_to = "PValue_Type", values_to = "PValue") %>%
  filter(!(Inheritance == "additive" & PValue_Type == "INHERENT_pvalue"))

# Perform ANOVA grouped by Inheritance_Model and PValue_Type
anova_results1 <- long_data1 %>%
  group_by(Inheritance) %>%
  do(tidy(aov(PValue ~ PValue_Type, data = .)))

# Perform Tukey's HSD test based on the corrected ANOVA setup
tukey_results1 <- long_data1 %>%
  group_by(Inheritance) %>%
  do({
    model <- aov(PValue ~ PValue_Type, data = .)
    tukey <- TukeyHSD(model)
    tidy(tukey)
  })

# Write to disk
write.csv(averages_filtered1, "Binary_AvgP_MAF01.csv", row.names = FALSE)
write.csv(anova_results1, "Binary_ANOVA_MAF01.csv", row.names = FALSE)
write.csv(tukey_results1, "Binary_Tukey_MAF01.csv", row.names = FALSE)

## Within Inheritance Model for MAF = 0.1 and PEN_DIFF = 0.1
Discrete_filtered2 <- Discrete_filtered1  %>% filter(PEN_DIFF == 0.1)

# Calculate average p-vals for each inheritance Model
averages_filtered2 <- Discrete_filtered2 %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_Original = mean(Original_pvalue, na.rm = TRUE),
    Avg_of_INHERENT = mean(INHERENT_pvalue, na.rm = TRUE),
    Avg_of_EDGE = mean(EDGE_pvalue, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_pvalue, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data2 <- pivot_longer(Discrete_filtered2, cols = c("Original_pvalue", "INHERENT_pvalue", "EDGE_pvalue", "PAGER_pvalue"),
                          names_to = "PValue_Type", values_to = "PValue") %>%
  filter(!(Inheritance == "additive" & PValue_Type == "INHERENT_pvalue"))

# Perform ANOVA grouped by Inheritance_Model and PValue_Type
anova_results2 <- long_data2 %>%
  group_by(Inheritance) %>%
  do(tidy(aov(PValue ~ PValue_Type, data = .)))

# Perform Tukey's HSD test based on the corrected ANOVA setup
tukey_results2 <- long_data2 %>%
  group_by(Inheritance) %>%
  do({
    model <- aov(PValue ~ PValue_Type, data = .)
    tukey <- TukeyHSD(model)
    tidy(tukey)
  })

# Write to disk
write.csv(averages_filtered2, "Binary_AvgP_MAF01_PD01.csv", row.names = FALSE)
write.csv(anova_results2, "Binary_ANOVA_MAF01_PD01.csv", row.names = FALSE)
write.csv(tukey_results2, "Binary_Tukey_MAF01_PD01.csv", row.names = FALSE)