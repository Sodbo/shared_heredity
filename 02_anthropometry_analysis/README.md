# Anthropometry analysis 
Anthropometry analysis is aimed to find shared heredity of four anthropometry traits: Body mass index (gwas_id=4049), Weight (gwas_id=4050), 
 Hip circumference (gwas_id=4058), Waist circumference (gwas_id=4179)
All traits uploaded to prod version of GWAS-Map

## Usage
 
The numeric prefix of scripts defines the order of using. Scripts without prefix number are used by other scripts and shouldn't be launched. Some of the scripts are used scripts from 00_core_function directory. Some scripts use GWAS-Map db functions and should be run from gwas-master (or other appropriate environment) 

## 00_start.sh
Starts the conveyer pipeline_for_calculation_of_matrices, and necessary for saving parametrs of running and to simplify running. The first argument is the path to the result directory. Names of result files is standard to automate passing its names as arguments of command line to following scripts.
Conveyer starts the following scripts:
01_pheno_corr.sh
02_convert_long_to_wide_form.R
03_gene_corr.sh
04_gene_corr_to_matrices.R
05_shared_heredity.R

