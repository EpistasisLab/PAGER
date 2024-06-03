# Load the required libraries
library(dplyr)
library(psych)
library(ggplot2)

setwd("/path/to/folder")

# Import SNPwise files
rep1 <- read.csv("/path/to/folder/replicate_1.csv")
rep2 <- read.csv("/path/to/folder/replicate_2.csv")
rep3 <- read.csv("/path/to/folder/replicate_3.csv")
rep4 <- read.csv("/path/to/folder/replicate_4.csv")
rep5 <- read.csv("/path/to/folder/replicate_5.csv")
rep6 <- read.csv("/path/to/folder/replicate_6.csv")
rep7 <- read.csv("/path/to/folder/replicate_7.csv")
rep8 <- read.csv("/path/to/folder/replicate_8.csv")
rep9 <- read.csv("/path/to/folder/replicate_9.csv")
rep10 <- read.csv("/path/to/folder/replicate_10.csv")

# Put all replicate dataframes in a list
rep_list <- list(rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep8, rep9, rep10)

# Remove the first column as it isn't part of the calculation
rep_list_nofirstcol <- lapply(rep_list, function(x) x[, -1])

# Average across dataframes in the list
average_rep <- Reduce("+", rep_list_nofirstcol) / length(rep_list_nofirstcol)

# Get back the SNP counts and reduce by 1
average_rep$SNPnum <- rep1$No.of.SNPs.
average_rep$SNPnum <- average_rep$SNPnum - 1

# Remove unwanted columns
average_rep$EDGE_TIME <- NULL
average_rep$PAGER_CPU_TIME <- NULL
average_rep$PAGER_GPU_TIME <- NULL

# Makee SNPnum first column
average_rep <- average_rep[, c(ncol(average_rep), 1:(ncol(average_rep)-1))]

# Add the PAGER computational time increase
average_rep$PAGER_CPU_INC <- average_rep$EDGE_TIME_WITHOUT_SNP1/average_rep$PAGER_CPU_TIME_WITHOUT_SNP1
average_rep$PAGER_GPU_INC <- average_rep$EDGE_TIME_WITHOUT_SNP1/average_rep$PAGER_GPU_TIME_WITHOUT_SNP1

# Write file to disk
write.csv(average_rep, "SNPwiseCompute.csv", row.names = FALSE)

## Discrete Plotting
# Fit a linear model for CPU
lm_CPU <- lm(PAGER_CPU_INC ~ SNPnum, data = average_rep)
lm_GPU <- lm(PAGER_GPU_INC ~ SNPnum, data = average_rep)

# Extract the coefficients
slope_CPU <- coef(lm_CPU)[["SNPnum"]]
intercept_CPU <- coef(lm_CPU)[["(Intercept)"]]
slope_GPU <- coef(lm_GPU)[["SNPnum"]]
intercept_GPU <- coef(lm_GPU)[["(Intercept)"]]

# Get logarithmic X breaks for X axis
x_min <- min(average_rep$SNPnum)
x_max <- max(average_rep$SNPnum)
min_pow <- floor(log10(x_min))
max_pow <- ceiling(log10(x_max))
breaks_at_powers_of_ten <- 10^seq(min_pow, max_pow, by = 1)

# Plot
ggplot(average_rep, aes(x = SNPnum)) +
  geom_point(aes(y = PAGER_CPU_INC), color = "black") +
  geom_point(aes(y = PAGER_GPU_INC), color = "black") +
  geom_smooth(aes(y = PAGER_CPU_INC), method = "lm", se = FALSE, color = "red") +
  geom_smooth(aes(y = PAGER_GPU_INC), method = "lm", se = FALSE, color = "blue") +
  theme_minimal() +
  labs(x = "Number of SNPs", y = "Speed Factor Increase") + 
  ggtitle("PAGER CPU/GPU Speed Increase over EDGE as SNPs Increase") +  
  scale_y_continuous(limits = c(3, 7), breaks = seq(3, 7, by = 1), expand = c(0, 0)) +
  scale_x_log10(breaks = breaks_at_powers_of_ten, labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_text(x = 1.9, y = 4.15, label = paste("CPU:","y = ", round(slope_CPU, 10), "x + ", round(intercept_CPU, 4)), hjust = 0, vjust = 0) +
  geom_text(x = 1.9, y = 5.3, label = paste("GPU:","y = ", round(slope_GPU, 10), "x + ", round(intercept_GPU, 4)), hjust = 0, vjust = 0) +
  theme(
    axis.title = element_text(size = 12, face = "bold"), 
    axis.text = element_text(size = 11), 
    plot.title = element_text(size = 14, face = "bold") 
  )

ggsave("SNPwise_ComputeTime_vs_EDGE.pdf", plot = last_plot(), device = "pdf", width = 10, height = 7, units = "in", useDingbats=FALSE)
