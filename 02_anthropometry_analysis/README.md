This folder contains scripts for anthropometry analysis aimed to find shared genetic background of five anthropometric traits: Body mass index (gwas_id=4049), Weight (gwas_id=4050), Hip circumference (gwas_id=4058), Waist circumference (gwas_id=4179) and Whole body fat mass (gwas_id=4054). All traits uploaded to prod version of GWAS-MAP database.

## Note
SH is an abbreviation of "shared heredity". Used in the code as the synonym of SGCT.

## Usage
The numeric prefix of scripts defines the order of using. Scripts without prefix number are used by other scripts and shouldn't be launched. Some of the scripts are used scripts from 00_core_function directory. Some scripts use GWAS-MAP dtabase functions and should be run from gwas-master (or other appropriate environment).
Numeration like 08_1b should be interpreted as following: 08 is order of using, 1b means variant 1 of analysis is for GIP1 (0 means SGCT and should be runned the first, 2 means maximal heretability or maxH) and b is the second step of validity checking.

## 00_reupload_gwases_to_test-db.sh
This script is optional and should be used if summary statistics file is not loaded to the shared_heredity container (GWAS-MAP database).

## 00_start.sh
Starts the conveyer pipeline_for_calculation_of_matrices, and necessary for saving parametrs of running and to simplify running. The first argument is the path to the result directory. Names of result files is standard to automate passing its names as arguments of command line to following scripts.
Conveyer starts the following scripts:
01_pheno_corr.sh
02_convert_long_to_wide_form.R
03_gene_corr.sh
04_gene_corr_to_matrices.R
05_shared_heredity.R

## 01_pheno_corr.sh
This script calculates phenotypic correlations between traits. The script is writted for GWAS-MAP database and it uses command line variables.

## 02_convert_long_to_wide_form.R
This script converts phenotypic correlations matrix from long to wide format. The script uses command line variables.

## 03_gene_corr.sh
This script calculates genetic correlations between traits. The script is writted for GWAS-MAP database and it uses command line variables.

## 04_gene_corr_to_matrices.R
This script calculates genetic correlations, covariances and se of genetic covariances between traits and builds matrices. The script uses command line variables.

## 05_alpha_coefficients.R
This script calculates alpha coeffitients and weights for each trait in linear combination. The script uses command line variables.

## 06_0_alpha_CI_estimation.R
This script estimates confident intervals for alpha coeffitients for SGCT.

## 06_1_alpha_CI_estimation_gip1.R
This script estimates confident intervals for alpha coeffitients for GIP1 (the GIP approach is described here doi: 10.1038/s42003-020-1051-9).

## 06_2_alpha_CI_estimation_maxH.R
This script estimates confident intervals for alpha coeffitients for MaxH (doi: 10.1159/000381641).

## 07_0_linear_combination.R
This script calculates summary statistics for SGCT.

## 07_1_linear_combination_gip1.R
This script calculates summary statistics GIP1.

## 07_2_linear_combination_maxH.R
This script calculates summary statistics for MaxH.

## 08_0_upload_sh_to_db.sh
This script uploads SGCT GWAS to the GWAS-MAP database.

## 08_1_upload_gip1_to_db.sh
This script uploads GIP1 GWAS to the GWAS-MAP database.

## 08_2_upload_maxH_to_db.sh
This script uploads MaxH GWAS to the GWAS-MAP database.

## 09_0b_validity_check.R
This script provides analytical estimation of heritability of SGCT and its genetic correlations with original traits using core functions.

## 09_1_validity_check_gip1.R
This is a script for internal use. It compares some earlier obtained results for GIP1 with the latest one.

## 09_1a_validity_check_gip1.sh
This script estimates GIP1 genetic correlations with original traits using LD Score regression implemented in GWAS-MAP database.

## 09_1b_validity_check_gip1.R
This script provides analytical estimation of GIP1 genetic correlations with original traits using core functions.

## 09_2_validity_check_maxH.R
This is a script for internal use. It compares some earlier obtained results for MaxH with the latest one.

## 09_2a_validity_check_maxH.sh
This script estimates MaxH genetic correlations with original traits using LD Score regression implemented in GWAS-MAP database.

## 09_2b_validity_check_maxH.R
This script provides analytical estimation of MaxH genetic correlations with original traits using core functions.

## 10_0_check_sh_intercept_for_gc_correction.sh
This script estimates heritability of SGCT and UGCTs using LD Score regression implemented in GWAS-MAP database. Intercept estimates are further used to decide whether the correction for genomic control is needed.

## 10_1_check_sh_intercept_for_gc_correction_gip1.sh
This script estimates heritability GIP1 using LD Score regression implemented in GWAS-MAP database. Intercept estimates are further used to decide whether the correction for genomic control is needed.

## 10_2_check_sh_intercept_for_gc_correction_maxH.sh
This script estimates heritability of the MaxH using LD Score regression implemented in GWAS-MAP database. Intercept estimate is further used to decide whether the correction for genomic control is needed.

## 11_0_p_val_gc_correction.R
This script performs a genomic control correction for original traits and SGCT.

## 11_1_p_val_gc_correction_gip1.R
This script performs a genomic control correction for GIP1.

## 11_2_p_val_gc_correction_maxH.R
This script performs a genomic control correction for MaxH.

## 13_linear_combination.r
This script adjusts original traits for SGCT.

## 14_upload_tr-sh_to_test.sh
This script uploads UGCTs summary statistics to the GWAS-MAP database.

## 15_start.sh
This script sources pipeline_for_calculation_of_matrices.sh script and sets the command line variables for further calculations (original traits, SGCT and UGCTs are mentioned).

## pipeline_for_calculation_of_matrices.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits, SGCT and UGCTs (see steps 01 - 05). The script uses command line variables.

## 16_validity_check.R
This script provides analytical estimation of heritability of traits adjusted for SGCT and their genetic correlations with SGCT using core functions.

## 16a_matrix_visualisation.R
This script creates a heatmap plot of genetic correlations of original traits, SGCT and UGCTs.

## 17_0_joint_clumping_and_enrichment.R
This script runs clumping for original traits and SGCT and finds shared hits.

## 17a_0_joint_clumping_and_enrichment.R
This script runs clumping for original traits, SGCT and UGCTs and finds shared hits.

## 18_0_joint_clumping_and_enrichment_gip1.R
This script runs clumping for original traits and GIP1 and finds shared hits.

## 19_0_joint_clumping_and_enrichment_maxH.R
This script runs clumping for original traits and MaxH and finds shared hits.

## 20_significant_loci.R
This script counts genome-wide significant loci for each original trait or SGCT.

## 21_significant_loci_sh.R
This script counts loci genome-wide significant both for particular original trait and for SGCT.

## DEPICT
This folder contains config files and script for DEPICT analysis for all anthropometric traits.

## test_upperW
This folder contains test script running the shared_heredity.R core function with upper W values.










