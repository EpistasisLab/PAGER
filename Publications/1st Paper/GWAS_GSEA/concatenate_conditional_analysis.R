# This is the script to concatenate all the CSV files for each encoding and each phenotype for the conditional analysis results

# Load necessary library
library(data.table)

# Set the working directory to where your CSV files are located
# Replace "/path/to/your/csv/files" with the actual path to your CSV files
mainDir <- "/common/compbiomed-dsn/projects/Philip_Freda_GWAS_20231024/analysis/conditional_analysis/pager_gwas_with_additive_GRMs_addsig/retrofat/"
setwd(mainDir)

# List all CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")

# Initialize an empty data table to store the concatenated data
all_data <- rbindlist(lapply(csv_files, fread), use.names = TRUE, fill = TRUE)

# Write the concatenated data to a new CSV file
output_file <- "retrofat_pager_additive_addsig_concatenated_data.csv"
fwrite(all_data, output_file)

cat("All CSV files have been concatenated into", output_file)
