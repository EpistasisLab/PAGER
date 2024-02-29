import os
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

# This function will take in one SNP and a Phenotype, will return the PAGER encoding values for the genotypes 0,1 and 2.
def pager_encoding_cpu(snp, phenotype):
    # create a dataframe with the single SNP and the phenotype
    snp_df = pd.DataFrame({'genotype': snp.astype(float), 'phenotype': phenotype.astype(float)})

    # Identify unique genotypic classes, excluding NAs
    unique_genotypes = snp_df['genotype'].dropna().unique()

    # Check if any of the expected genotypic classes (0, 1, and 2) are missing
    missing_genotypes = set([0, 1, 2]) - set(unique_genotypes)

    # Handle missing genotypic classes
    if missing_genotypes:
        # Choose the first available minimum genotype as the new anchor
        anchor = min(unique_genotypes)
    else:
        anchor = 0

    # Recode the genotypes based on their mean phenotype values, only for genotypes 0, 1, and 2
    geno_aggregations = snp_df[(snp_df['genotype'].isin([0, 1, 2])) & (~snp_df['genotype'].isna())].groupby('genotype').agg(
        mean_phenotype=('phenotype', 'mean')
    )

    # add the genotype values to the geno_aggregations dataframe
    geno_aggregations['genotype'] = geno_aggregations.index

    # Calculate the relative distances from the mean phenotype of the anchor genotype
    anchor_mean = geno_aggregations[geno_aggregations['genotype'] == anchor]['mean_phenotype'].values[0]
    geno_aggregations['rel_dist'] = (geno_aggregations['mean_phenotype'] - anchor_mean) 

    # Use Min-Max normalization on rel_dist
    scaler = MinMaxScaler()
    geno_aggregations['normalized_rel_dist'] = scaler.fit_transform(geno_aggregations['rel_dist'].values.reshape(-1, 1))

    # return geno_aggregations which is a dataframe with the PAGER encoding values, it will have three columns, genotype, mean_phenotype and rel_dist
    return geno_aggregations

def process_file(input_file_path):

    data_to_encode = pd.read_csv(input_file_path)

    phenotype_file = pd.read_csv("/home/ghosha/common/bams_edge_pager/encoding_for_gwas/phenotype_bmi.csv")

    # Create a new dataframe to store the encoded values
    pager_encoded_data = pd.DataFrame()

    # Iterate over each SNP column in data_to_encode
    for col in data_to_encode.columns:
        # Check if all values in the column are NA
        if data_to_encode[col].isna().all():
            # If all values are NA, add the column as it is to the new dataframe
            pager_encoded_data[col] = data_to_encode[col]
        else:
            # Otherwise, apply pager_encoding_cpu to the SNP column and phenotype
            snp_encoding = pager_encoding_cpu(data_to_encode[col], phenotype_file['bmi'])

            # Use map to replace the original SNP values with the new encoding values
            pager_encoded_data[col] = data_to_encode[col].map(snp_encoding.set_index('genotype')['normalized_rel_dist'])

    # Extract the file name without extension
    file_name = os.path.splitext(os.path.basename(input_file_path))[0]

    # Save the results to a CSV file with the edge_encoded prefix
    output_folder = "/home/ghosha/common/bams_edge_pager/encoding_for_gwas/pager_encoded_permuted"
    output_file_path = f"{output_folder}/pager_encoded_{file_name}.csv"
    pager_encoded_data.to_csv(output_file_path, index=False)

def main():

    # Parse command-line argument for the input file
    import sys
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <input_file>")
        sys.exit(1)

    input_file_path = sys.argv[1]

    process_file(input_file_path)

if __name__ == "__main__":
    main()
