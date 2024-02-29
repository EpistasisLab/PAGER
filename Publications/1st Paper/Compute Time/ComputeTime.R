# Load the required libraries
library(dplyr)
library(psych)
library(ggplot2)

setwd("/Users/fredap/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/all simulation tests/single_snp_discrete/Figures/ComputeTime")

# Import discrete files
add_disc <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/additive_discrete.csv")
dom_disc <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/dominant_discrete.csv")
rec_disc <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/recessive_discrete.csv")
subadd_disc <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/subadditive_discrete.csv")
supadd_disc <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/superadditive_discrete.csv")
het_disc <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/heterosis_discrete.csv")
und_disc <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/underdominant_discrete.csv")
ovd_disc <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/overdominant_discrete.csv")

# Import continuous files
add_cont <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/additive_continuous.csv")
dom_cont <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/dominant_continuous.csv")
rec_cont <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/recessive_continuous.csv")
subadd_cont <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/subadditive_continuous.csv")
supadd_cont <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/superadditive_continuous.csv")
het_cont <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/heterosis_continuous.csv")
und_cont <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/underdominant_continuous.csv")
ovd_cont <- read.csv("~/Library/CloudStorage/Box-Box/CedarsSinai/AutoQTL/EDGE/Results/For manuscript/cpu vs gpu/Redone_with_formula_change/overdominant_continuous.csv")

## Discrete
# Modify dataframes and add an Inheritance column:
# Define the columns you want to extract
columns_to_extract <- c("Num_Samples", "Case_Control_Ratio", "MAFA", "PEN_DIFF", "EDGE_TIME", "PAGER_CPU_TIME", "PAGER_GPU_TIME")

# Modify each data frame to contain only the info we want
add_disc <- add_disc[(columns_to_extract)]
dom_disc <- dom_disc[(columns_to_extract)]
rec_disc <- rec_disc[(columns_to_extract)]
subadd_disc <- subadd_disc[(columns_to_extract)]
supadd_disc <- supadd_disc[(columns_to_extract)]
het_disc <- het_disc[(columns_to_extract)]
und_disc <- und_disc[(columns_to_extract)]
ovd_disc <- ovd_disc[(columns_to_extract)]

# Add the inheritance column
add_disc$Inheritance <- "Additive"
dom_disc$Inheritance <- "Dominant"
rec_disc$Inheritance <- "Recessive"
subadd_disc$Inheritance <- "Subadditive"
supadd_disc$Inheritance <- "Superadditive"
het_disc$Inheritance <- "Heterosis"
und_disc$Inheritance <- "Underdominant"
ovd_disc$Inheritance <- "Overdominant"

# Create a list of the data frames
data_frames_disc <- list(add_disc, dom_disc, rec_disc, subadd_disc, supadd_disc, het_disc, und_disc, ovd_disc)

# Combine all modified data frames into one
combined_disc <- do.call(rbind, data_frames_disc)

# Extract just times and sample size because we are averaging across all inheritance patterns and other settings
times_disc <- combined_disc[,c(1,5,6,7)]

# Outlier control to remove initialization of GPU
# Calculate mean and standard deviation for the specified column
mean_value <- mean(times_disc$PAGER_GPU_TIME)
std_dev <- sd(times_disc$PAGER_GPU_TIME)

# Define a Z-score threshold (e.g., 2 standard deviations)
z_score_threshold <- 3

# Filter the dataframe to keep only rows within the Z-score threshold
times_disc <- times_disc[abs((times_disc$PAGER_GPU_TIME - mean_value) / std_dev) <= z_score_threshold, ]

# Group the data by Num_Samples and calculate the means and standard deviations for PAGER and EDGE times
result_disc <- times_disc %>%
  group_by(Num_Samples) %>%
  summarise(
    EDGE_Mean = mean(EDGE_TIME),
    EDGE_SD = sd(EDGE_TIME) / sqrt(n()),
    PAGER_CPU_Mean = mean(PAGER_CPU_TIME),
    PAGER_CPU_SD = sd(PAGER_CPU_TIME) / sqrt(n()),
    PAGER_GPU_Mean = mean(PAGER_GPU_TIME),
    PAGER_GPU_SD = sd(PAGER_GPU_TIME) / sqrt(n())
  )

# Add the method column
result_disc$Method <- "Discrete"

# Add the PAGER computational time increases for CPU and GPU
result_disc$PAGER_CPU_INC <- result_disc$EDGE_Mean/result_disc$PAGER_CPU_Mean
result_disc$PAGER_GPU_INC <- result_disc$EDGE_Mean/result_disc$PAGER_GPU_Mean

## Add in the continuous data
# Modify dataframes and add an Inheritance column:
# Define the columns you want to extract
columns_to_extract <- c("Num_Samples", "MAFA", "PEN_DIFF", "EDGE_TIME", "PAGER_CPU_TIME", "PAGER_GPU_TIME")

# Modify each data frame to contain only the info we want
add_cont <- add_cont[(columns_to_extract)]
dom_cont <- dom_cont[(columns_to_extract)]
rec_cont <- rec_cont[(columns_to_extract)]
subadd_cont <- subadd_cont[(columns_to_extract)]
supadd_cont <- supadd_cont[(columns_to_extract)]
het_cont <- het_cont[(columns_to_extract)]
und_cont <- und_cont[(columns_to_extract)]
ovd_cont <- ovd_cont[(columns_to_extract)]

# Add the inheritance column
add_cont$Inheritance <- "Additive"
dom_cont$Inheritance <- "Dominant"
rec_cont$Inheritance <- "Recessive"
subadd_cont$Inheritance <- "Subadditive"
supadd_cont$Inheritance <- "Superadditive"
het_cont$Inheritance <- "Heterosis"
und_cont$Inheritance <- "Underdominant"
ovd_cont$Inheritance <- "Overdominant"

# Create a list of the data frames
data_frames_cont <- list(add_cont, dom_cont, rec_cont, subadd_cont, supadd_cont, het_cont, und_cont, ovd_cont)

# Combine all modified data frames into one
combined_cont <- do.call(rbind, data_frames_cont)

# Extract just times and sample size because we are averaging across all inheritance patterns and other settings
times_cont <- combined_cont[,c(1,4,5,6)]

# Outlier control to remove initialization of GPU
# Calculate mean and standard deviation for the specified column
mean_value <- mean(times_cont$PAGER_GPU_TIME)
std_dev <- sd(times_cont$PAGER_GPU_TIME)

# Define a Z-score threshold (e.g., 3 standard deviations)
z_score_threshold <- 3

# Filter the dataframe to keep only rows within the Z-score threshold
times_cont <- times_cont[abs((times_cont$PAGER_GPU_TIME - mean_value) / std_dev) <= z_score_threshold, ]

# Group the data by Num_Samples and calculate the means and standard deviations for PAGER and EDGE times
result_cont <- times_cont %>%
  group_by(Num_Samples) %>%
  summarise(
    EDGE_Mean = mean(EDGE_TIME),
    EDGE_SD = sd(EDGE_TIME) / sqrt(n()),
    PAGER_CPU_Mean = mean(PAGER_CPU_TIME),
    PAGER_CPU_SD = sd(PAGER_CPU_TIME) / sqrt(n()),
    PAGER_GPU_Mean = mean(PAGER_GPU_TIME),
    PAGER_GPU_SD = sd(PAGER_GPU_TIME) / sqrt(n())
  )

# Add the method column
result_cont$Method <- "Continuous"

# Add the PAGER computational time increase
result_cont$PAGER_CPU_INC <- result_cont$EDGE_Mean/result_cont$PAGER_CPU_Mean
result_cont$PAGER_GPU_INC <- result_cont$EDGE_Mean/result_cont$PAGER_GPU_Mean

## Discrete Plotting
# Fit a linear model for CPU
lm_CPU <- lm(PAGER_CPU_INC ~ Num_Samples, data = result_disc)
lm_GPU <- lm(PAGER_GPU_INC ~ Num_Samples, data = result_disc)

# Extract the coefficients
slope_CPU <- coef(lm_CPU)[["Num_Samples"]]
intercept_CPU <- coef(lm_CPU)[["(Intercept)"]]
slope_GPU <- coef(lm_GPU)[["Num_Samples"]]
intercept_GPU <- coef(lm_GPU)[["(Intercept)"]]

ggplot(result_disc, aes(x = Num_Samples)) +
  geom_point(aes(y = PAGER_CPU_INC), color = "black") +
  geom_point(aes(y = PAGER_GPU_INC), color = "black") +
  geom_smooth(aes(y = PAGER_CPU_INC), method = "lm", se = FALSE, color = "red") +
  geom_smooth(aes(y = PAGER_GPU_INC), method = "lm", se = FALSE, color = "blue") +
  theme_minimal() +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +  
  labs(x = "Sample Size", y = "Speed Factor Increase") + 
  ggtitle("PAGER CPU/GPU Speed Increase over EDGE (Discrete Phenotype)") +  
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 50000), breaks = seq(0, 50000, by = 5000), expand = c(0.03, 0)) +
  geom_text(x = max(as.numeric(result_disc$Num_Samples)), y = max(result_disc$PAGER_CPU_INC), label = paste("CPU:","y = ", round(slope_CPU, 5), "x + ", round(intercept_CPU, 4)), hjust = 1.5, vjust = 14) +
  geom_text(x = max(as.numeric(result_disc$Num_Samples)), y = max(result_disc$PAGER_GPU_INC), label = paste("GPU:","y = ", round(slope_GPU, 5), "x + ", round(intercept_GPU, 4)), hjust = 1.5, vjust = 4) +
  theme(
    axis.title = element_text(size = 12, face = "bold"), 
    axis.text = element_text(size = 11), 
    plot.title = element_text(size = 14, face = "bold") 
  )

ggsave("Discrete_ComputeTime_vs_EDGE.pdf", plot = last_plot(), device = "pdf", width = 10, height = 7, units = "in", useDingbats=FALSE)

## Continuous Plotting
# Fit a linear model for CPU
lm_CPU_cont <- lm(PAGER_CPU_INC ~ Num_Samples, data = result_cont)
lm_GPU_cont <- lm(PAGER_GPU_INC ~ Num_Samples, data = result_cont)

# Extract the coefficients
slope_CPU_cont <- coef(lm_CPU_cont)[["Num_Samples"]]
intercept_CPU_cont <- coef(lm_CPU_cont)[["(Intercept)"]]
slope_GPU_cont <- coef(lm_GPU_cont)[["Num_Samples"]]
intercept_GPU_cont <- coef(lm_GPU_cont)[["(Intercept)"]]

ggplot(result_cont, aes(x = Num_Samples)) +
  geom_point(aes(y = PAGER_CPU_INC), color = "black") +
  geom_point(aes(y = PAGER_GPU_INC), color = "black") +
  geom_smooth(aes(y = PAGER_CPU_INC), method = "lm", se = FALSE, color = "red") +
  geom_smooth(aes(y = PAGER_GPU_INC), method = "lm", se = FALSE, color = "blue") +
  theme_minimal() +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +  
  labs(x = "Sample Size", y = "Speed Factor Increase") + 
  ggtitle("PAGER CPU/GPU Speed Increase over EDGE (Continuous Phenotype)") +  
  scale_y_continuous(limits = c(0, 41), breaks = seq(0, 41, by = 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 50000), breaks = seq(0, 50000, by = 5000), expand = c(0.03, 0)) +
  geom_text(x = max(as.numeric(result_cont$Num_Samples)), y = max(result_cont$PAGER_CPU_INC), label = paste("CPU:","y = ", round(slope_CPU_cont, 5), "x + ", round(intercept_CPU_cont, 4)), hjust = 1.5, vjust = 12) +
  geom_text(x = max(as.numeric(result_cont$Num_Samples)), y = max(result_cont$PAGER_GPU_INC), label = paste("GPU:","y = ", round(slope_GPU_cont, 5), "x + ", round(intercept_GPU_cont, 4)), hjust = 1.5, vjust = 4) +
  theme(
    axis.title = element_text(size = 12, face = "bold"), 
    axis.text = element_text(size = 11), 
    plot.title = element_text(size = 14, face = "bold") 
  )

ggsave("Cont_ComputeTime_vs_EDGE.pdf", plot = last_plot(), device = "pdf", width = 10, height = 7, units = "in", useDingbats=FALSE)
