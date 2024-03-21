# This script is for permuting the bimbam files for additive, dominant and recessive

# load libraries
import pandas as pd
import numpy as np
import os
import sys

# Give the complete path of the directory which contains the bimbam files you want to permute
path = sys.argv[1]

def shuffle_ids(df):
    fixed_cols = df.iloc[:, :3]  # Keep the first three columns unchanged
    genotype_cols = df.iloc[:, 3:]  # Columns to be shuffled
    # Ensure NA values are shuffled correctly
    shuffled_genotype_cols = genotype_cols.apply(lambda x: np.random.permutation(x.to_numpy()), axis=0)
    return pd.concat([fixed_cols, shuffled_genotype_cols], axis=1)

chromosome_folder = path

chromosome_files = os.listdir(chromosome_folder)

# All shuffled files are saved in a folder - shuffled_data_0001 is all chromosome bimbam files for the first shuffle and so on
for i in range(1000):
    shuffle_folder = f'shuffled_data_{i:04d}'
    os.makedirs(shuffle_folder, exist_ok=True)
    np.random.seed(i)

    for file in chromosome_files:
        file_path = os.path.join(chromosome_folder, file)
        # Ensure NA values are read correctly
        genotype_data = pd.read_csv(file_path, header=None, delim_whitespace=True, na_values='NA')
        shuffled_data = shuffle_ids(genotype_data)
        # Add a check to ensure the number of columns remains consistent
        assert shuffled_data.shape[1] == genotype_data.shape[1], "Column mismatch after shuffling"
        shuffled_file_path = os.path.join(shuffle_folder, f'shuffled_{file}')
        # Ensure the file is written with the same format
        shuffled_data.to_csv(shuffled_file_path, sep=' ', index=False, header=False, na_rep='NA')
