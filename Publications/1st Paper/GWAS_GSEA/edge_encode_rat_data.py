from pandas_genomics import sim, scalars
from pandas_genomics.arrays import GenotypeArray, GenotypeDtype

import concurrent.futures
import pandas as pd
import numpy as np
from enum import Enum
import time
import random
from sklearn.preprocessing import MinMaxScaler
import os
import clarite
import time

import statsmodels.api as sm
import warnings
warnings.filterwarnings("ignore")

# Name conversation for parameters
def conv_u(i):
    switcher = {
        "REC": "RECESSIVE",
        "SUB": "SUB_ADDITIVE",
        "ADD": "ADDITIVE",
        "SUP": "SUPER_ADDITIVE",
        "DOM": "DOMINANT",
        "HET": "HET",
        "NUL": "ADDITIVE",
        "OVD": "OVERDOMINANT",
        "UND": "UNDERDOMINANT",
    }
    return switcher.get(i, "Invalid Group")


def conv_l(i):
    switcher = {
        "REC": "Recessive",
        "SUB": "Sub-Additive",
        "ADD": "Additive",
        "SUP": "Super_Additive",
        "DOM": "Dominant",
        "HET": "Heterozygous",
        "NUL": "NULL",
        "OVD": "OVERDOMINANT",
        "UND": "UNDERDOMINANT",
    }
    return switcher.get(i, "Invalid Group")

variant1 = scalars.Variant("1", 1, id="rs1", ref="A", alt=["C"])
variant2 = scalars.Variant("1", 2, id="rs2", ref="G", alt=["T"])

def get_snp1_gt_array(gt_table_idxs):
        """Assemble a GenotypeArray for SNP1 directly from genotype table indices"""
        # Check for missing values (NaN) in gt_table_idxs
        has_missing_values = np.isnan(gt_table_idxs).any()
    
        if has_missing_values:
            print("gt_table_idxs contains missing values (NaN).")

        print(gt_table_idxs)
        snp1 = scalars.Variant("1", 1, id="rs1", ref="A", alt=["C"])
        dtype = GenotypeDtype(snp1)
        gt_table_data = (
            ((0, 0), 255),
            ((0, 1), 255),
            ((1, 1), 255),
            ((0, 0), 255),
            ((0, 1), 255),
            ((1, 1), 255),
            ((0, 0), 255),
            ((0, 1), 255),
            ((1, 1), 255)
        )
        print("i values:", [i for i in gt_table_idxs])
        data = np.array(
            [gt_table_data[i] for i in gt_table_idxs], dtype=dtype._record_type
        )
        return GenotypeArray(values=data, dtype=dtype)


def sim_bootstrap(discrete, ab1,ab2,case_control_ratio,num_samples,PEN_DIFF,MAFA,train_seed):
    n_controls = int(num_samples*(1-case_control_ratio))
    n_cases = num_samples - n_controls
    PEN_BASE = (1-PEN_DIFF)/2
    MAFB = MAFA

    sim_name = "ab1_"+str(ab1)+"_casecontrolrat_"+str(case_control_ratio) + "_pendiff_" + str(PEN_DIFF) + "_maf_" + str(MAFA)
    # print(sim_name)

    ALL_RESULTS_ENCODING = pd.DataFrame()
    ALL_RESULTS_EDGE_ALPHA = pd.DataFrame()

    # print(ab1, ab2)
    ab1u = conv_u(ab1)
    ab2u = conv_u(ab2)
    ab1l = conv_l(ab1)
    # BLOCK 1
    # Main Effect for SNP1 without interaction
    # Training data
    train_main = sim.BAMS.from_model(
        eff1=getattr(sim.SNPEffectEncodings, ab1u),
        eff2=getattr(sim.SNPEffectEncodings, ab2u),
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=train_seed)

    # BLOCK 2
    if discrete:
        train_me = train_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB)
        train_me = train_me.sort_values(by="Outcome",ascending=False)
        
    else:
        train_me = train_main.generate_quantitative(n_samples=num_samples, maf1=MAFA, maf2=MAFB)
        train_me = train_me.sort_values(by="Outcome",ascending=False)

    train_me = train_me.reset_index(drop=True)
    
    np.random.seed(train_seed)
    Age_1 = pd.DataFrame(np.random.poisson(lam=31, size=(int(0.507*n_controls))).astype(int))
    Age_2 = pd.DataFrame(np.random.poisson(lam=54.5, size=(int(0.316*n_controls))).astype(int))
    Age_3 = pd.DataFrame(np.random.poisson(lam=82.5, size=(n_controls-int(0.507*n_controls)-int(0.316*n_controls))).astype(int))
    Age_4 = pd.DataFrame(np.random.poisson(lam=31, size=(int(0.144*n_cases))).astype(int))
    Age_5 = pd.DataFrame(np.random.poisson(lam=54.5, size=(int(0.435*n_cases))).astype(int))
    Age_6 = pd.DataFrame(np.random.poisson(lam=82.5, size=(n_cases-int(0.144*n_cases)-int(0.435*n_cases))).astype(int))
    Age = pd.concat([Age_1,Age_2,Age_3,Age_4,Age_5,Age_6], axis=0, ignore_index=True, sort=False)
    BMI_1 = pd.DataFrame(np.random.normal(loc=26.5, scale=1.0, size=(n_controls)).astype(int))
    BMI_2 = pd.DataFrame(np.random.normal(loc=35, scale=7.5, size=(n_cases)).astype(int))
    BMI = pd.concat([BMI_1,BMI_2], axis=0, ignore_index=True, sort=False)
    Sex_1 = pd.DataFrame(np.random.binomial(n=1, p=0.481, size=(n_controls)).astype(int))
    Sex_2 = pd.DataFrame(np.random.binomial(n=1, p=0.525, size=(n_cases)).astype(int))
    Sex = pd.concat([Sex_1,Sex_2], axis=0, ignore_index=True, sort=False)
    Smoking_1 = pd.DataFrame(np.random.binomial(n=1, p=0.125, size=(n_controls)).astype(int))
    Smoking_2 = pd.DataFrame(np.random.binomial(n=1, p=0.64, size=(n_cases)).astype(int))
    Smoking = pd.concat([Smoking_1,Smoking_2], axis=0, ignore_index=True, sort=False) 
    train_COV = pd.concat([Age,BMI,Sex,Smoking], axis=1, ignore_index=True, sort=False)
    train_cov = pd.concat([train_me, train_COV.reindex(train_me.index)], axis=1)
    train_cov.columns=['Outcome','SNP1','SNP2','Age', 'BMI', 'Sex','Smoking']
    train_cov = train_cov[['Outcome','SNP1','Age', 'BMI', 'Sex','Smoking']]

    return train_cov, train_me


#############################################################################
# starting the whole pipeline of analysis

def process_file(input_file_path):

    # code to encode a single data file
    seed = random.randint(0, 100000)
    train_seed = seed

    phenotype_file = pd.read_csv("/home/ghosha/common/bams_edge_pager/encoding_for_gwas/phenotype_bmi.csv")

    # # Parse command-line argument for the input file
    # import sys
    # if len(sys.argv) != 2:
    #     print("Usage: python your_script.py <input_file>")
    #     sys.exit(1)

    # input_file_path = sys.argv[1]

    # Your existing code for reading data and processing
    data_to_encode = pd.read_csv(input_file_path)

    print("Data to encode: ", data_to_encode)

    edge_encoded_data = pd.DataFrame()
    #edge_encoded_data['V1'] = data_to_encode['']

    for column in data_to_encode.columns[0:]:
        print("Column: ", column)
        edge_data = data_to_encode[[column]].copy()
        edge_data['bmi'] = phenotype_file['bmi']

        # Check if all values in the column are NaN
        if pd.isna(edge_data[column]).all():
            # If all values are NaN, copy the column to edge_encoded_data
            edge_encoded_data[column] = edge_data[column]
            continue  # Skip the rest of the loop for this column


        num_samples = data_to_encode[column].notna().sum()
        if num_samples != data_to_encode.shape[0]:
            missing_values = 'True'
        else:
            missing_values = 'False'
        work_group = (
            ["REC", "ADD"],
            ["SUB", "ADD"],
            ["ADD", "ADD"],
            ["SUP", "ADD"],
            ["DOM", "ADD"],
            ["HET", "ADD"],
            ["OVD", "ADD"],
            ["UND", "SUB"],
            ["NUL", "ADD"])

        case_control_ratio = 1/4 # [1/4, 1/2]
        pen_diff = 0.4 #[0.05, 0.1, 0.175, 0.25, 0.33, 0.4]:
        MAFA = 0.4 # [0.05, 0.1, 0.2, 0.3, 0.4]:

        ab1, ab2 = work_group[2]

        discrete = True

        num_snps = 1

        train_cov, train_me = sim_bootstrap(discrete, ab1,ab2,case_control_ratio,num_samples,pen_diff,MAFA,train_seed)

        features = edge_data[column]
        target = edge_data['bmi']

        train_me['Outcome'] = target

        feature_array = features.values
        #print("i values:", [i for i in feature_array])

        # Check for missing values (NaNs) and print their indices
        missing_value_indices = np.where(np.isnan(feature_array))[0]
        print("Missing value indices:", missing_value_indices)          
            
        # Removing NA values
        # Create a mask for non-NaN values
        non_nan_mask = ~np.isnan(feature_array)
        # Filter feature_array to keep only non-NaN values
        feature_array = feature_array[non_nan_mask]

        # remove the corresponding rows having NAs from target as well
        target = target[non_nan_mask]
        print("Shape of target: ", target.shape)

        # Check for missing values (NaNs) and print their indices after removing NAs
        missing_value_indices = np.where(np.isnan(feature_array))[0]
        print("Missing value indices after removing NAs:", missing_value_indices)     

        # Convert non-NaN values to integers, keep NaN values as NaN
        feature_array = feature_array.astype(int)
        print("Feature array values after converting to int:", [i for i in feature_array])
            
        geno_data = get_snp1_gt_array(feature_array)
        train_me["SNP1"] = geno_data
        edge_weights_me_t = train_me.genomics.calculate_edge_encoding_values(
        data=train_me, outcome_variable="Outcome", covariates=[])
        alpha_value = edge_weights_me_t.loc[0, 'Alpha Value']

        edge_data[column] = edge_data[column].replace({1: alpha_value, 2: 1})

        print("EDGE DATA: ", edge_data)
        edge_encoded_data[column] = edge_data[column]
            
    # # add the bmi at the end of the file
    # edge_encoded_data['bmi'] = phenotype_file['bmi']
    # Save the results to a CSV file
    # Extract the file name without extension
    file_name = os.path.splitext(os.path.basename(input_file_path))[0]

    # Save the results to a CSV file with the edge_encoded prefix
    output_folder = "/home/ghosha/common/bams_edge_pager/encoding_for_gwas/edge_encoded_permuted"
    output_file_path = f"{output_folder}/edge_encoded_{file_name}.csv"
    edge_encoded_data.to_csv(output_file_path, index=False)


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