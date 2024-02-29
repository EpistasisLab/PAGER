# This is the script to get the lowest pvalues after running GWAS on shuffled encoded data

import pandas as pd
import os

# Define the base path and the encodings
base_path = "/common/compbiomed-dsn/projects/Philip_Freda_GWAS_20231024/"
#encodings = ["additive_gwas", "dominant_with_additive_GRMs", "recessive_with_additive_GRMs"]
#encodings = ["edge_gwas"]
encodings = ["pager_gwas"]

# Iterate over each encoding
for encoding in encodings:
    # Define the path for the current encoding's shuffled data
    #shuffled_data_path = os.path.join(base_path, encoding, "bmiwtail/permutations/slurm_scripts_for_gwas/output/")
    #shuffled_data_path = os.path.join(base_path, encoding, "edge_gwas_with_additive_GRMs/bmiwtail/second_permutations/slurm_scripts_for_gwas/output/")
    shuffled_data_path = os.path.join(base_path, encoding, "pager_gwas_with_additive_GRMs/bmiwtail/second_permutations/slurm_scripts_for_gwas/output/")
    
    # Initialize a list to store the results for this encoding
    results = []

    # Iterate over each shuffled_data folder, from 0000 to 0999
    for folder_number in range(100):  # This includes 0000 to 0999
        folder_name = f"shuffled_data_{folder_number:04d}"
        folder_path = os.path.join(shuffled_data_path, folder_name)

        # Initialize a variable to store the lowest p_score for the current folder
        lowest_p_score = None

        # Iterate over each chromosome file
        for chromosome in range(1, 21):
            file_path = os.path.join(folder_path, f"chr{chromosome}.gwas.out.assoc.txt")

            # Read the file and extract the lowest p_score value
            try:
                df = pd.read_csv(file_path, delim_whitespace=True)
                min_p_score = df['p_score'].min()

                # Update the lowest_p_score if it's the lowest so far
                if lowest_p_score is None or min_p_score < lowest_p_score:
                    lowest_p_score = min_p_score
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

        # Append the result for the current folder
        results.append([folder_number, lowest_p_score])

    # Convert the results to a DataFrame and save as a CSV file
    results_df = pd.DataFrame(results, columns=['Folder', 'Lowest_p_score'])
    results_csv_path = os.path.join(shuffled_data_path, f'lowest_p_score_values_{encoding}_second_permutations.csv')
    results_df.to_csv(results_csv_path, index=False)

    print(f"CSV file created for {encoding}: {results_csv_path}")

# End of script
print("All CSV files created.")
