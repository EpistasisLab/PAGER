# This is the script to run GWAS on the shuffled encoded data

import os

# Define the base path and the encodings
base_path = "/common/compbiomed-dsn/projects/Philip_Freda_GWAS_20231024/"
encodings = ["pager_gwas/pager_gwas_with_additive_GRMs/bmiwtail/second_permutations"]

# Loop through each encoding
for encoding in encodings:
    encoding_path = os.path.join(base_path, encoding)
    additive_path = os.path.join(encoding_path, "additive_GRMs")
    snp_annotation_path = os.path.join(encoding_path, "snp_annotation_files")
    slurm_scripts_path = os.path.join(encoding_path, "slurm_scripts_for_gwas")
    shuffled_bimbam_base_folder = os.path.join(encoding_path, "shuffled_csv_files")  # Base path for shuffled BIMBAM files
    output_path = os.path.join(slurm_scripts_path, "output")

    # Create directories if they don't exist
    if not os.path.exists(slurm_scripts_path):
        os.makedirs(slurm_scripts_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Loop through each folder (numbered from 0000 to 0099)
    for folder_number in range(0, 100):
        convert_folder_name = f"convert_{folder_number:04d}"
        shuffled_bimbam_folder = os.path.join(shuffled_bimbam_base_folder, convert_folder_name)
        folder_name = f"shuffled_data_{folder_number:04d}"
        folder_path = os.path.join(output_path, folder_name)

        # Create a results sub-folder for the current folder
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        # Initialize SLURM script content for the entire folder
        script_content = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -p defq
#SBATCH --job-name=gwas_{convert_folder_name}
#SBATCH -t 12:00:00
#SBATCH --mem=50GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=priyanka.bhandary@cshs.org
#SBATCH -e {slurm_scripts_path}/{folder_name}.err
#SBATCH -o {slurm_scripts_path}/{folder_name}.out

module load singularity-apptainer/1.1.6

"""

        # Loop through each chromosome (1 to 20)
        for chromosome in range(1, 21):
            genome_file = os.path.join(shuffled_bimbam_folder, f"chr{chromosome}_pager_encoded_shuffled_data_{folder_number:04d}.bimbam")
            phenotype_file = os.path.join(encoding_path, "pheno_final.txt")
            kinship_matrix_file = os.path.join(additive_path, f"allExcept.chr{chromosome}_additive.GRM.cXX.txt")
            snp_annotation_file = os.path.join(snp_annotation_path, f"chr{chromosome}_snp_annotation_file")

            # Append the command for this chromosome to the script content
            script_content += f"""
singularity exec -B {encoding_path}:{encoding_path} docker://pgcbioinfo/gemma gemma \
-g {genome_file} \
-p {phenotype_file} \
-n 1 \
-k {kinship_matrix_file} \
-a {snp_annotation_file} \
-lmm 4 \
-o {folder_name}/chr{chromosome}.gwas.out
"""

        # Write the combined script to a file in the slurm_scripts directory
        script_filename = f"{slurm_scripts_path}/run_gwas_{convert_folder_name}.slurm"
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

        print(f"Created combined script: {script_filename}")

# End of script
print("All combined SLURM scripts created.")
