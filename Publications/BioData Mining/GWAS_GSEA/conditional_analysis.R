# This is the script to conduct conditional analysis

#load libraries
library(data.table)

# Set working directory
mainDir <- "/common/compbiomed-dsn/projects/Philip_Freda_GWAS_20231024/analysis/conditional_analysis/"
setwd(mainDir)

# Setting command line arguments
args <- commandArgs(trailingOnly = TRUE)

trait <- args[1] # bmiwtail/bodyweight/retrofat
chr <- args[2] # chr1-chr20
gwas_threshold<-as.numeric(args[3]) # Threshold from permutation significance testing results
encoding<-args[4] # additive/dominant_additive/recessive_additive/pager/edge
output_path<-args[5] # output path
kinship_file<-args[6] # GRM file
phenotype_file<-args[7] # phenotype file
snp_annotation_file<-args[8] # SNP annotation file
bimbam_file_path<-args[9] #bimbam files

print(gwas_threshold)
# convert GWAS to numeric
#gwas_threshold<-as.numeric(gwas_threshold)

# Read in files
covariate_file_base_path=paste0(output_path,"temp/") # Read in covariate file for each chromosome for each trait called in arguments - we are doing this per chromosome
gwas_summary_base_path=paste0(output_path,"gwas_summary_files/") # GWAS files obtained

# Initialize variable and data structures
cov_counter <- 0 # Set the covariate counter to 0 
loop_run <- TRUE # Set loop to TRUE
top_snps_info <- list() # Initilaize a list for the top SNPs

# Begin the main loop for conditional analysis
while (loop_run) {
  # Define the GWAS summary statistics file path for the current iteration. For the first iteration, read in the gwas summary stats file from the initial run. If any other iteration, read in the previous run which will always be saved in the "ouput" folder in the working directory
  gwas_summary_path <- ifelse(cov_counter == 0,
                              paste0(gwas_summary_base_path, chr, ".gwas.out.assoc.txt"),
                              paste0("/common/compbiomed-dsn/projects/Philip_Freda_GWAS_20231024/analysis/conditional_analysis/output/", "output_gwas_iteration_", cov_counter-1, "_", encoding, "_", trait, ".assoc.txt"))
  
  # Load GWAS summary statistics
  # Load GWAS summary statistics
  gwas_data <- fread(gwas_summary_path)
  
  gwas_data <- gwas_data[p_score != 0]
  # Identify the top SNP based on -log10 P-value
  top_snp <- gwas_data[which.min(p_score)]
  print(top_snp)
  
  # Resolve duplicates if any and select the SNP with the highest absolute effect size
  if (sum(gwas_data$p_score == top_snp$p_score, na.rm = TRUE) > 1) {
    duplicates <- gwas_data[gwas_data$p_score == top_snp$p_score]
    top_snp <- duplicates[which.max(abs(duplicates$beta))]
  }
  
  # Update top SNP information
  top_snps_info[[length(top_snps_info) + 1]] <- top_snp
  print(top_snps_info)

  # Read in the bimbam file
  bimbam_file_path_chr <- paste0(bimbam_file_path, chr, "_", encoding, ".bimbam")
  
  system(paste("grep", top_snp$rs, bimbam_file_path_chr, "> extracted_snp_genotypes.txt"), intern = FALSE)
  snp_genotype_data <- scan("extracted_snp_genotypes.txt", what = "", sep = "\n")
  
  # Assuming the genotype data is numeric and starts from the 4th column in the BIMBAM file
  top_snp_genotypes <- strsplit(snp_genotype_data, " ")[[1]][-(1:3)]
  
  # Prepare the covariate file path
  covariate_file_path <- sprintf("%scovariates_%s_iter_%d.txt", covariate_file_base_path, chr, cov_counter)
  
  # Prepare the covariate file
  covariate_file_path <- sprintf("%scovariates_%s_iter_%d.txt", covariate_file_base_path, chr, cov_counter)
  # For the first iteration, create the initial covariate file
  if (cov_counter == 0) {
    # Prepare initial covariate structure. It might include only the intercept or the intercept + first significant SNP.
    # Here, let's start with just the intercept and add the first SNP's genotype data.
    accumulated_covariates <- data.frame(Intercept = rep(1, 3166))  # Adjust this number to match your dataset size
    new_snp_column_name <- "SNP_1"  # Naming the first SNP column
    accumulated_covariates[[new_snp_column_name]] <- top_snp_genotypes  # Add the first SNP genotype data
  } else {
    # For subsequent iterations, add new SNP genotype data to accumulated_covariates
    new_snp_column_name <- paste("SNP", cov_counter + 1, sep = "_")
    accumulated_covariates[[new_snp_column_name]] <- top_snp_genotypes
  }
  
  # Write the accumulated_covariates to the covariate file, ensuring column names are included
  write.table(accumulated_covariates, covariate_file_path, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  
  # Check if this is the first run and if the top SNP's p-value is greater than the threshold
  if (cov_counter == 0 && top_snp$p_score > gwas_threshold) {
    cat("The p-value of the top SNP exceeds the GWAS threshold on the first iteration. Stopping the analysis.\n")
    loop_run <- FALSE  # This will stop the while loop
    break  # In case the loop_run flag isn't checked immediately, this ensures the script stops here
  }
  
  # Generate a SLURM script for this iteration - this will change according to your HPC specifications
  slurm_script_path <- sprintf("%s/gemma_slurm_iter_%d.sh", output_path, cov_counter)
  output_file_prefix <- sprintf("output_gwas_iteration_%d_%s_%s", cov_counter, encoding, trait)
  gemma_command <- sprintf("singularity exec -B %s:%s docker://pgcbioinfo/gemma gemma -g %s -p %s -n 1 -k %s -a %s -lmm 4 -c %s -o %s",
                           output_path, output_path, bimbam_file_path_chr, phenotype_file, kinship_file, snp_annotation_file,
                           covariate_file_path, output_file_prefix)
  print(gemma_command)
  singularity_load <- sprintf("module load singularity-apptainer/1.1.6")
  slurm_script_content <- sprintf(
    "#!/bin/bash\n#SBATCH --job-name=gemma_analysis_%d\n#SBATCH --output=%s/gemma_output_%%j.txt\n#SBATCH --error=%s/gemma_error_%%j.txt\n#SBATCH --time=02:00:00\n#SBATCH --mem=8G\n#SBATCH --cpus-per-task=2\n%s\n%s",
    cov_counter, output_path, output_path, singularity_load, gemma_command)
  
  # Write the SLURM script to a file
  writeLines(slurm_script_content, slurm_script_path)
  
  # Submit the SLURM script and capture job ID
  job_submission_output <- system(paste("sbatch", slurm_script_path), intern = TRUE)
  job_id <- gsub("Submitted batch job ", "", job_submission_output)
  
  # Assuming a simple repeat loop for checking job status
  repeat {
    Sys.sleep(60)  # Check every minute
    running_jobs <- system("squeue -u $USER", intern = TRUE)  # List your jobs in the queue
    if (!any(grepl(job_id, running_jobs))) {
      break  # Exit the loop if the job ID is no longer listed
    }
  }
  
  # After ensuring the job has completed, proceed to load new GWAS results
  new_gwas_summary_path <- sprintf("/common/compbiomed-dsn/projects/Philip_Freda_GWAS_20231024/analysis/conditional_analysis/output/output_gwas_iteration_%d_%s_%s.assoc.txt", cov_counter, encoding, trait)
  # Begin the main loop for conditional analysis
  if (file.exists(new_gwas_summary_path)) {
    new_gwas_data <- fread(new_gwas_summary_path)
    # Check for additional significant SNPs
    if (all(new_gwas_data$p_score > gwas_threshold & new_gwas_data$p_score != 0, na.rm = TRUE)) {
      print(gwas_threshold)
      print(min(new_gwas_data$p_score))
      save(top_snps_info, file = paste0(output_path, "/", trait, "_", encoding, "_", chr, "_conditional_analysis.RData"))
      write.csv(rbindlist(top_snps_info, use.names = TRUE, fill = TRUE), paste0(output_path, "/", trait, "_", encoding, "_", chr, "_conditional_analysis.csv"), row.names = FALSE)
      loop_run <- FALSE
    } else {
      print(min(new_gwas_data$p_score))
      print(gwas_threshold)
      cov_counter <- cov_counter + 1
    }
  } else {
    warning("GWAS summary file does not exist after job completion:", new_gwas_summary_path)
    loop_run <- FALSE  # Consider breaking the loop if the expected output file is missing
  }
}



