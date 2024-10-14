### R
### PJF
### 3/14/24

### This script will show average p-values for a host of certain combinations within inheritance models

## Continuous
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
additive <- additive[, c(1,2,3,4,5,6,7,8,9,14,10,11,12,13)]

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
combined_cont <- bind_rows(dataframes, .id = "Inheritance")

# Multiply EDGE pvalue by 2 to control for multiple testing
combined_cont$EDGE_pvalue <- combined_cont$EDGE_pvalue*2

## All Data
# Calculate average p-vals for each inheritance Model
averages_filtered <- combined_cont %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_Original = mean(Original_pvalue, na.rm = TRUE),
    Avg_of_INHERENT = mean(INHERENT_pvalue, na.rm = TRUE),
    Avg_of_EDGE = mean(EDGE_pvalue, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_pvalue, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data <- pivot_longer(combined_cont, cols = c("Original_pvalue", "INHERENT_pvalue", "EDGE_pvalue", "PAGER_pvalue"),
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
setwd("path/to/folder")
write.csv(averages_filtered, "Cont_AvgP_All.csv", row.names = FALSE)
write.csv(anova_results, "Cont_ANOVA_All.csv", row.names = FALSE)
write.csv(tukey_results, "Cont_Tukey_All.csv", row.names = FALSE)

## Within Inheritance Model for MAF = 0.1
Cont_filtered1 <- combined_cont %>% filter(MAFA == 0.1)

# Calculate average p-vals for each inheritance Model
averages_filtered1 <- Cont_filtered1 %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_Original = mean(Original_pvalue, na.rm = TRUE),
    Avg_of_INHERENT = mean(INHERENT_pvalue, na.rm = TRUE),
    Avg_of_EDGE = mean(EDGE_pvalue, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_pvalue, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data1 <- pivot_longer(Cont_filtered1, cols = c("Original_pvalue", "INHERENT_pvalue", "EDGE_pvalue", "PAGER_pvalue"),
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
setwd("/path/to/folder")
write.csv(averages_filtered1, "Cont_AvgP_MAF01.csv", row.names = FALSE)
write.csv(anova_results1, "Cont_ANOVA_MAF01.csv", row.names = FALSE)
write.csv(tukey_results1, "Cont_Tukey_MAF01.csv", row.names = FALSE)

## Within Inheritance Model for MAF = 0.1 and PEN_DIFF = 0.1
Cont_filtered2 <- Cont_filtered1  %>% filter(PEN_DIFF == 0.1)

# Calculate average p-vals for each inheritance Model
averages_filtered2 <- Cont_filtered2 %>%
  group_by(Inheritance) %>%
  summarise(
    Avg_of_Original = mean(Original_pvalue, na.rm = TRUE),
    Avg_of_INHERENT = mean(INHERENT_pvalue, na.rm = TRUE),
    Avg_of_EDGE = mean(EDGE_pvalue, na.rm = TRUE),
    Avg_of_PAGER = mean(PAGER_pvalue, na.rm = TRUE)
  )

# Convert data from wide to long format, excluding INHERENT_pvalue for additive model
long_data2 <- pivot_longer(Cont_filtered2, cols = c("Original_pvalue", "INHERENT_pvalue", "EDGE_pvalue", "PAGER_pvalue"),
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
setwd("/path/to/folder")
write.csv(averages_filtered2, "Cont_AvgP_MAF01_PD01.csv", row.names = FALSE)
write.csv(anova_results2, "Cont_ANOVA_MAF01_PD01.csv", row.names = FALSE)
write.csv(tukey_results2, "Cont_Tukey_MAF01_PD01.csv", row.names = FALSE)