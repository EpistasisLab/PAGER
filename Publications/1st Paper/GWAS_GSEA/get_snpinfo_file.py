# This is the script to get the SNP annotation file for conducting GWAS

import os, sys

input_filename = sys.argv[1] # This is the bimbam file for the specific chromosome
snpinfo_filename = sys.argv[2] # This is the output SNP annotation file

# Read and write files
input_file = open(input_filename, 'r')
snpinfo_file = open(snpinfo_filename, 'w')

# From the bimbam file, save the chromosome and base pair information to a dictionary
snpid_to_info = {}
for line in input_file:
    snpid = line.strip().split()
    snpid = snpid[0]
    #print(snpid)
    chromosome = snpid.split(".")[0]
    bp_info = snpid.split(".")[1]
    snpid_to_info[snpid] = [chromosome, bp_info]
    
#print(snpid_to_info)

# Write the SNP, base pair information and chromosome to the output file
for snpid in snpid_to_info:
    snpinfo_file.write(snpid+", "+snpid_to_info[snpid][1]+", "+snpid_to_info[snpid][0])
    snpinfo_file.write("\n")