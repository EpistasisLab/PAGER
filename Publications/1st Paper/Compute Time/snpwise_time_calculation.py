"""This script is used to read csv files of varying number of SNPs (10, 100, 1000, 10000, 1000000) with a continuous phenotype and then 
finding the time required to calculate the EDGE encoding value, PAGER encoding value using CPU and GPU. """

# all the import statements
import os
import cupy as cp
import pandas as pd
import numpy as np
from pandas_genomics import sim, scalars
from sklearn.preprocessing import MinMaxScaler
import statsmodels.api as sm
import random
import clarite
import time
from pandas_genomics.arrays import GenotypeArray, GenotypeDtype
import math

############################### BAMS DATA SIMULATION ##############################

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


def simulate_genotype(discrete, ab1,ab2,case_control_ratio,num_samples,PEN_DIFF,MAFA,train_seed):
    n_controls = int(num_samples*(1-case_control_ratio))
    n_cases = num_samples - n_controls
    PEN_BASE = (1-PEN_DIFF)/2
    MAFB = MAFA

    sim_name = "ab1_"+str(ab1)+"_casecontrolrat_"+str(case_control_ratio) + "_pendiff_" + str(PEN_DIFF) + "_maf_" + str(MAFA)

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
    else:
        train_me = train_main.generate_quantitative(n_samples=num_samples, maf1=MAFA, maf2=MAFB)

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

    #starting the timer to calculate EDGE encoding time
    edge_start_time = time.time()  
    edge_weights_me_t = train_me.genomics.calculate_edge_encoding_values(
        data=train_cov, outcome_variable="Outcome",covariates=['Age', 'BMI', 'Sex','Smoking'])
    # ending the timer
    edge_end_time = time.time()
    elapsed_time = edge_end_time - edge_start_time
    edge_weights_me = edge_weights_me_t.copy()
    edge_weights_me.insert(loc=0, column="BioAct", value=ab1l)
    edge_weights_me.insert(loc=0, column="TrainSeed", value=train_seed)

    #change column type to string
    train_cov['SNP1'] = train_cov['SNP1'].astype(str)
    # encode the SNP1 to 0,1,2 from A/A, A/C, C/C
    train_cov['SNP1'] = train_cov['SNP1'].replace("A/A", 0)
    train_cov['SNP1'] = train_cov['SNP1'].replace('A/C', 1)
    train_cov['SNP1'] = train_cov['SNP1'].replace('C/C', 2)
    
    #rearrange the columns put SNP1 and Outcome at the end
    train_cov = train_cov[['Age', 'BMI', 'Sex','Smoking','SNP1','Outcome']]
    train_cov.columns = ['Age', 'BMI', 'Sex','Smoking','SNP','Phenotype']
    if discrete:
        train_cov['Phenotype'] = train_cov['Phenotype'].replace("Case", 1)
        train_cov['Phenotype'] = train_cov['Phenotype'].replace("Control", 0)

    return train_cov, edge_weights_me_t, elapsed_time

###################################################################################################

########### PAGER CPU FUNCTION ###############

# This function will take in one SNP and a Phenotype, will return the PAGER encoding values for the genotypes 0,1 and 2.
def pager_encoding_cpu(snp, phenotype):
    # create a dataframe with the single SNP and the phenotype
    snp_df = pd.DataFrame({'genotype': snp.astype(float), 'phenotype': phenotype.astype(float)})

    # calculate the phenotypic mean per genotype
    geno_aggregations = snp_df.groupby('genotype').agg(
        mean_phenotype = ('phenotype', 'mean')
    )

    # add the genotype values to the geno_aggregations dataframe
    geno_aggregations['genotype'] = geno_aggregations.index

    # use the PAGER formula
    AA_mean = geno_aggregations.loc[geno_aggregations['genotype'].idxmin(), 'mean_phenotype']
    geno_aggregations['rel_dist'] = (geno_aggregations['mean_phenotype'] - (AA_mean))

    # Use Min-Max normalization on rel_dist
    scaler = MinMaxScaler()
    geno_aggregations['normalized_rel_dist'] = scaler.fit_transform(geno_aggregations['rel_dist'].values.reshape(-1, 1))

    # return geno_aggregations which is a dataframe with the PAGER encoding values, it will have three columns, genotype, mean_phenotype and rel_dist
    return geno_aggregations


###################################################################################

def pager_encoding_gpu(snp, phenotype):
    snp = cp.asarray(snp, dtype=cp.float32)
    phenotype = cp.asarray(phenotype, dtype=cp.float32)

    # Get unique types in the 'types' array
    unique_snps = cp.unique(snp)

    # Calculate the mean of each type
    means = [cp.mean(phenotype[snp == t]) for t in unique_snps]

    means_arrays = [arr.get() for arr in means]

    geno_aggregations = pd.DataFrame({'genotype': cp.asnumpy(unique_snps), 'mean_phenotype': means_arrays})

    geno_aggregations['genotype'] = geno_aggregations.index

    # Use the PAGER formula
    AA_mean = geno_aggregations.loc[geno_aggregations['genotype'].idxmin(), 'mean_phenotype']
    
    # changed formula without std
    geno_aggregations['rel_dist'] = (geno_aggregations['mean_phenotype'] - AA_mean)

    # Use Min-Max normalization on rel_dist
    # Find the minimum and maximum values in the NumPy array
    min_value = np.min(geno_aggregations['rel_dist'])
    max_value = np.max(geno_aggregations['rel_dist'])

    # Scale and shift the array to fit between -1 and 1
    scaled_array = -1 + 2 * (geno_aggregations['rel_dist'] - min_value) / (max_value - min_value)

    geno_aggregations['normalized_rel_dist'] = scaled_array
    return geno_aggregations


def get_snp1_gt_array(gt_table_idxs):
    """Assemble a GenotypeArray for SNP1 directly from genotype table indices"""
    cleaned_idxs = []

    for idx in gt_table_idxs:
        if not math.isnan(idx):
            cleaned_idxs.append(int(idx))
        else:
            cleaned_idxs.append(np.nan)

    print(type(cleaned_idxs[0]))
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
        ((1, 1), 255),
    )

    # Create an array with a placeholder value, e.g., -1
    data = np.full(len(cleaned_idxs), -1, dtype=dtype._record_type)
    
    # Replace elements with valid values and keep NaN elements as NaN
    for i, idx in enumerate(cleaned_idxs):
        if not math.isnan(idx):
            data[i] = gt_table_data[int(idx)]

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


#############################

# Read the real data and get the edge, pager_cpu and pager_gpu time
seed = random.randint(0, 100000)
train_seed = seed
# Create an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['No of SNPs.', 'EDGE_TIME', 'EDGE_TIME_WITHOUT_SNP1', 'PAGER_CPU_TIME', 'PAGER_CPU_TIME_WITHOUT_SNP1', 'PAGER_GPU_TIME', 'PAGER_GPU_TIME_WITHOUT_SNP1'])


data_folder = "/path/to/folder"

# Iterate through all the files in the folder to get the edge value, then pager encoding
for filename in os.listdir(data_folder):
    if filename.endswith(".csv"):
        print("Name of the file: ", filename)
        edge_total_time = 0
        pager_cpu_total_time = 0
        pager_gpu_total_time = 0
        edge_total_time_without_snp1 = 0
        pager_cpu_total_time_without_snp1 = 0
        pager_gpu_total_time_without_snp1 = 0

        # Extract "No of SNPs" and "Set no." from the file name
        split_parts = filename.split("_")
        no_of_snps = split_parts[-1]

        # Read the original data
        data = pd.read_csv(os.path.join(data_folder, filename))
        num_snps = 0

        for column in data.columns[:-1]:
            original_data = data[[column]].copy()
            original_data['bmi'] = data['bmi']
            edge_data = data[[column]].copy()
            edge_data['bmi'] = data['bmi']
            pager_cpu_data = data[[column]].copy()
            pager_cpu_data['bmi'] = data['bmi']
            pager_gpu_data = data[[column]].copy()
            pager_gpu_data['bmi'] = data['bmi']

            num_snps = num_snps + 1
        
            num_samples = data[column].notna().sum()
            if num_samples != data.shape[0]:
                missing_values = 'True'
            else:
                missing_values = 'False'
            # Define your parameter settings
            work_group = (
                ["REC", "ADD"],
                ["SUB", "ADD"],
                ["ADD", "ADD"],
                ["SUP", "ADD"],
                ["DOM", "ADD"],
                ["HET", "ADD"],
                ["OVD", "ADD"],
                ["UND", "ADD"],
                ["NUL", "ADD"]
                )

            # using a default assignment as no simulation is required
            case_control_ratio = 1/2 
            pen_diff = 0.4 
            MAFA = 0.4 
            
            ab1, ab2 = work_group[2]

            discrete = True


            train_cov, train_me = sim_bootstrap(discrete, ab1,ab2,case_control_ratio,num_samples,pen_diff,MAFA,train_seed)

            features = edge_data[column]
            target = edge_data['bmi']

            train_me['Outcome'] = target

            feature_array = features.values

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
            edge_start_time = time.time()
            edge_weights_me_t = train_me.genomics.calculate_edge_encoding_values(
                    data=train_me, outcome_variable="Outcome", covariates=[])
            edge_end_time = time.time()
            edge_time = edge_end_time - edge_start_time
            edge_total_time = edge_total_time + edge_time
            if num_snps > 1:
                edge_total_time_without_snp1 = edge_total_time_without_snp1 + edge_time

            print("EDGE TIME per SNP: ", edge_time)

            # Perform PAGER CPU encoding for the data
            pager_cpu_start_time = time.time()
            pager_encoding_result = pager_encoding_cpu(pager_cpu_data[column], pager_cpu_data['bmi'])
            pager_cpu_end_time = time.time()
            pager_cpu_time = pager_cpu_end_time - pager_cpu_start_time
            pager_cpu_total_time = pager_cpu_total_time + pager_cpu_time
            if num_snps > 1:
                pager_cpu_total_time_without_snp1 = pager_cpu_total_time_without_snp1 + pager_cpu_time
            print("PAGER CPU TIME per SNP: ", pager_cpu_time)

            # Perform PAGER GPU encoding for the data
            pager_gpu_start_time = time.time()
            pager_encoding_result = pager_encoding_gpu(pager_gpu_data[column], pager_gpu_data['bmi'])
            pager_gpu_end_time = time.time()
            pager_gpu_time = pager_gpu_end_time - pager_gpu_start_time
            pager_gpu_total_time = pager_gpu_total_time + pager_gpu_time
            if num_snps > 1:
                pager_gpu_total_time_without_snp1 = pager_gpu_total_time_without_snp1 + pager_gpu_time
            print("PAGER GPU TIME per SNP: ", pager_gpu_time)

    # Append the results to the DataFrame
    result_df = result_df.append({
        'No of SNPs.': no_of_snps,
        'EDGE_TIME': edge_total_time,
        'EDGE_TIME_WITHOUT_SNP1': edge_total_time_without_snp1,
        'PAGER_CPU_TIME': pager_cpu_total_time,
        'PAGER_CPU_TIME_WITHOUT_SNP1': pager_cpu_total_time_without_snp1,
        'PAGER_GPU_TIME': pager_gpu_total_time,
        'PAGER_GPU_TIME_WITHOUT_SNP1': pager_gpu_total_time_without_snp1
    }, ignore_index=True)

# Save the results to a CSV file
result_df.to_csv("/path/to/save/replicate_1.csv", index=False)
