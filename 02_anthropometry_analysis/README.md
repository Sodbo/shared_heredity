# Anthropometry analysis 
Anthropometry analysis is aimed to find shared heredity of five anthropometry traits: Body mass index (gwas_id=4049), Weight (gwas_id=4050), 
 Hip circumference (gwas_id=4058), Waist circumference (gwas_id=4179)
 Whole body fat mass (gwas_id=4054)
All traits uploaded to prod version of GWAS-Map

## Usage
 
The numeric prefix of scripts defines the order of using. Scripts without prefix number are used by other scripts and shouldn't be launched. Some of the scripts are used scripts from 00_core_function directory. Some scripts use GWAS-Map db functions and should be run from gwas-master (or other appropriate environment)
Numeration like 08_1b should be interpreted as following: 08 is order of using, 1b means variant 1 of analysis is for GIP1 (0 means shared heredity and should be runned the first, 2 means maximal heretability or maxH) and b is the second step of validity checking.

## 00_reupload_gwases_to_test-db.sh
is optional and should be used if summary statistics file is not loaded to shared heredity container

## 00_start.sh
Starts the conveyer pipeline_for_calculation_of_matrices, and necessary for saving parametrs of running and to simplify running. The first argument is the path to the result directory. Names of result files is standard to automate passing its names as arguments of command line to following scripts.
Conveyer starts the following scripts:
01_pheno_corr.sh
02_convert_long_to_wide_form.R
03_gene_corr.sh
04_gene_corr_to_matrices.R
05_shared_heredity.R

