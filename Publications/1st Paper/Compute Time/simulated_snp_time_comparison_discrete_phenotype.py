"""This script is used to compare the time required to calculate EDGE encoding values, PAGER encoding values using a CPU and GPU for 
SNPs with additive encoding and discrete phenotype."""

# all the import statements
import cupy as cp
import pandas as pd
import numpy as np
from pandas_genomics import sim, scalars
from sklearn.preprocessing import MinMaxScaler
import statsmodels.api as sm
import random
import time

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
        # train_me = train_me.sort_values(by="Outcome",ascending=False)
    else:
        train_me = train_main.generate_quantitative(n_samples=num_samples, maf1=MAFA, maf2=MAFB)
        # train_me = train_me.sort_values(by="Outcome",ascending=False)

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

    #starting the timer to find time required for calculate EDGE encoding value
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
    AA_mean = geno_aggregations.loc[0, 'mean_phenotype']
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
    AA_mean = geno_aggregations.loc[0, 'mean_phenotype']
    geno_aggregations['rel_dist'] = (geno_aggregations['mean_phenotype'] - AA_mean) 


    # Use Min-Max normalization on rel_dist
    # Find the minimum and maximum values in the NumPy array
    min_value = np.min(geno_aggregations['rel_dist'])
    max_value = np.max(geno_aggregations['rel_dist'])

    # Scale and shift the array to fit between -1 and 1
    scaled_array = -1 + 2 * (geno_aggregations['rel_dist'] - min_value) / (max_value - min_value)

    geno_aggregations['normalized_rel_dist'] = scaled_array
    return geno_aggregations

#############################

# Simulate data and then derive the edge time, pager cpu time and pager gpu time

# Define your parameter settings
work_group = (
    ["REC", "ADD"], # recessive encoding
    ["SUB", "ADD"], # subadditive encoding
    ["ADD", "ADD"], # additive encoding
    ["SUP", "ADD"], # superadditive encoding
    ["DOM", "ADD"], # dominant encoding
    ["HET", "ADD"], # heterosis encoding
    ["OVD", "ADD"], # overdominant encoding
    ["UND", "ADD"], # underdominant encoding
    ["NUL", "ADD"] # null encoding
)

discrete = True # type of phenotype
num_samples = [2000, 5000, 10000, 15000, 20000, 25000, 50000] # number of samples
PEN_DIFF_values = [0.1, 0.25, 0.33, 0.4] # penetrance difference
MAFA_values = [0.1, 0.2, 0.3, 0.4, 0.5] # minor allele frequency
case_control_ratio_values = [0.25, 0.5] # case control ratio

# Define the number of iterations
num_iterations = 1000 # number of replicates

# Create an empty DataFrame to store the results
results_df = pd.DataFrame(columns=["Iteration", "Num_Samples", "Case_Control_Ratio", "MAFA", "PEN_DIFF",  "EDGE_TIME", "PAGER_CPU_TIME", "PAGER_GPU_TIME"])

# Loop through the parameter settings and run the simulation
for iteration in range(num_iterations):
        for case_control_ratio in case_control_ratio_values:
            for num_samples_setting in num_samples:
                for PEN_DIFF_setting in PEN_DIFF_values:
                    for MAFA_setting in MAFA_values:
                        # Generate a random train seed for each iteration
                        train_seed = random.randint(1, 1000000)
                        ab1, ab2 = work_group[2] # decides the inheritence pattern of the SNPs (0 - recessive, 1 - subadditive, 2 - additive, 3 - superadditive, 4 - dominant, 5 - heterosis, 6 - overdominant, 7 - underdominant, 8 - null )
                        # Call your simulate_genotype function with the current settings
                        train_cov, edge_value_dataframe, edge_elapsed_time = simulate_genotype(
                            discrete,
                            ab1,
                            ab2,
                            case_control_ratio,
                            num_samples_setting,
                            PEN_DIFF_setting,
                            MAFA_setting,
                            train_seed
                        )
                        print("EDGE encoding time: ", edge_elapsed_time)

                        # a copy of the original data to be used for PAGER encoding
                        train_cov_copy = train_cov.copy()
                        
                        # PAGER CPU TIME CALCULATIOn
                        pager_start_time = time.time()
                        pager_encoding_result = pager_encoding_cpu(train_cov_copy['SNP'], train_cov_copy['Phenotype'])
                        # end timer for PAGER
                        pager_end_time = time.time()
                        pager_cpu_elapsed_time = pager_end_time - pager_start_time
                        print("PAGER CPU encoding time: ", pager_cpu_elapsed_time)

                        # PAGER GPU CALCULATION
                        pager_start_time = time.time()
                        pager_encoding_result = pager_encoding_gpu(train_cov_copy['SNP'], train_cov_copy['Phenotype'])
                        # end timer for PAGER
                        pager_end_time = time.time()
                        pager_gpu_elapsed_time = pager_end_time - pager_start_time
                        print("PAGER GPU encoding time: ", pager_gpu_elapsed_time)
                        
                        # Append the results to the DataFrame
                        results_df = results_df.append({
                            "Iteration": iteration,
                            "Num_Samples": num_samples_setting,
                            "Case_Control_Ratio": case_control_ratio,
                            "MAFA": MAFA_setting,
                            "PEN_DIFF": PEN_DIFF_setting,
                            "EDGE_TIME": edge_elapsed_time,
                            "PAGER_CPU_TIME": pager_cpu_elapsed_time,
                            "PAGER_GPU_TIME": pager_gpu_elapsed_time
                        }, ignore_index=True)

# Save the results DataFrame to a CSV file
results_df.to_csv("/home/ghosha/common/bams_edge_pager/cpu_vs_gpu/new_formula_additive_discrete.csv", index=False)

# Print the total execution time for calculate_edge_encoding_values
total_edge_encoding_time = sum(results_df['EDGE_TIME'])

# Print the total execution time for pager cpu
total_pager_cpu_encoding_time = sum(results_df['PAGER_CPU_TIME'])

# Print the total execution time for pager gpu
total_pager_gpu_encoding_time = sum(results_df['PAGER_GPU_TIME'])

print("Total EDGE encoding time: ", total_edge_encoding_time)
print("Total PAGER CPU encoding time: ", total_pager_cpu_encoding_time)
print("Total PAGER GPU encoding time: ", total_pager_gpu_encoding_time)
