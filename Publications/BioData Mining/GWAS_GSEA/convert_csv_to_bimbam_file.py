# This script is to convert the encodings in the CSV files to bimbam file format

import os
import sys
import csv

csv_filename = sys.argv[1] # Read in the CSV files with encoding
original_bimbam_filename = sys.argv[2] # The original bimbam files from Apurva that has the allele information Note: this script is run per chromosome so please use the specific bimbam file for the chromosome you are interested in
bimbam_filename = sys.argv[3] # The bimbam file that you want to create Note: This should be the same chromosome as the original bimbam file

# Open the files for reading and writing
csv_file = open(csv_filename, 'r')
bimbam_file = open(bimbam_filename, 'w')

# These function were written in case the SNPs in the original bimbam file and the CSV files have different formats - SNP:XXX or SNP.XXX
def preprocess_snp_id(snp_id):
    # Replace ":" with "." in SNP IDs
    return snp_id.replace(":", ".")

def postprocess_snp_id(snp_id):
    return snp_id.replace(".", ":")

# Step 1: Read the first file and store SNP data in a dictionary
snp_data = {}
with open(original_bimbam_filename, "r") as original_bimbam_file:
    for line in original_bimbam_file:
        parts = line.strip("\n").split() 
        snp_id = parts[0]
        allele1 = parts[1]
        allele2 = parts[2]
        snp_data[preprocess_snp_id(snp_id)] = [allele1, allele2] # Create a dictionary with the SNP as key and alleles as values

# Step 2: Read the second file to get column names
data_dict = {}
with open(csv_filename, "r") as csv_file:
    csv_reader = csv.reader(csv_file)

    header = next(csv_reader)

    # Initialize the dictionary with keys from the header
    for column_name in header:
        data_dict[column_name.strip()] = []

    # Read the rest of the data and populate the dictionary
    for row in csv_reader:
        for i, column_name in enumerate(header):
            # Replace blank spaces with "NA"
            if row[i].strip() == "":
                row[i] = "NA"
            data_dict[column_name].append(row[i])

# Step 3-4: Match SNP IDs, extract alleles, and row values
result_dict = {}

for snp_id, alleles in snp_data.items():
    if snp_id in data_dict:
        result_dict[snp_id] = {
            'Allele1': alleles[0],
            'Allele2': alleles[1],
            'RowValues': data_dict[snp_id]
        }

# Print or write the result_dict to the bimbam file
with open(bimbam_filename, 'w') as bimbam_file:
    for snp_id, data in result_dict.items():
        row_values = ' '.join(data['RowValues'])
        bimbam_file.write(f"{snp_id} {data['Allele1']} {data['Allele2']} {row_values}\n")