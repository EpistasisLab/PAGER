# This is the script to get the bimbam files from shuffled encoded data for PAGER and EDGE

import os

# Define the base path and the script path
base_path = "/common/compbiomed-dsn/projects/Philip_Freda_GWAS_20231024/pager_gwas/pager_gwas_with_additive_GRMs/bmiwtail/second_permutations/"
script_path = os.path.join(base_path, "scripts/convert_csv_to_bimbam_file.py")
original_bimbam_folder = os.path.join(base_path, "originial_bimbam_files")
output_bimbam_folder = os.path.join(base_path, "shuffled_csv_files")
slurm_scripts_folder = os.path.join(base_path, "slurm_scripts")

# Ensure the slurm_scripts directory exists
os.makedirs(slurm_scripts_folder, exist_ok=True)

# SLURM script template
slurm_script_template = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -p defq
#SBATCH --job-name=bimbam_{folder_name}
#SBATCH -t 12:00:00
#SBATCH --mem=50GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=priyanka.bhandary@cshs.org
#SBATCH --output={output_file}
#SBATCH --error={error_file}

"""

# Iterate over each shuffled data file
for file_number in range(100):  # from 0000 to 0099
    file_name = f"pager_encoded_shuffled_data_{file_number:04d}.csv"
    shuffled_file_path = os.path.join(output_bimbam_folder, file_name)
    folder_name = f"convert_{file_number:04d}"

    # Create a unique output folder for each script
    script_output_folder = os.path.join(output_bimbam_folder, folder_name)
    os.makedirs(script_output_folder, exist_ok=True)

    # Define the output and error file paths
    output_file = os.path.join(script_output_folder, f"{folder_name}.out")
    error_file = os.path.join(script_output_folder, f"{folder_name}.err")

    # Initialize SLURM script content
    slurm_script_content = slurm_script_template.format(
        folder_name=folder_name,
        output_file=output_file,
        error_file=error_file
    )

    # Create a command for each chromosome
    for chromosome in range(1, 21):
        original_bimbam_file = os.path.join(original_bimbam_folder, f"chr{chromosome}.round2_impute2_3473.bimbam")
        output_bimbam_file = os.path.join(script_output_folder, f"chr{chromosome}_pager_encoded_shuffled_data_{file_number:04d}.bimbam")

        command = f"python {script_path} {shuffled_file_path} {original_bimbam_file} {output_bimbam_file}\n"
        slurm_script_content += command

    # Write the SLURM script to a file
    slurm_script_filename = os.path.join(slurm_scripts_folder, f"run_convert_{file_number}.slurm")
    with open(slurm_script_filename, 'w') as script_file:
        script_file.write(slurm_script_content)

    print(f"Created SLURM script: {slurm_script_filename}")

# End of script
print("All SLURM scripts created.")
