This folder contains all code extracting shared heredity from psychometric traits.

## 00a_start.sh
This script sources 00b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (only original traits are mentioned). 

## 00b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits and estimation of alpha coeffitients for linear combination (see steps 01 - 05). The script uses command line variables.

## 01_pheno_corr.sh
This script calculates phenotypic correlations between traits. The script is writted for GWAS-Map database and it uses command line variables.

## 02_convert_long_to_wide_form.R
This script converts phenotypic correlations matrix from long to wide format. The script uses command line variables.

## 03_gene_corr.sh
This script calculates genetic correlations between traits. The script is writted for GWAS-Map database and it uses command line variables.

## 04_gene_corr_to_matrices.R
This script calculates genetic correlations, covariances and se of genetic covariances between traits and builds matrices. The script uses command line variables.

## 05_alpha_coefficients.R
This script calculates alpha coeffitients and weights for each trait in linear combination. The script uses command line variables.

## 06_alpha_CI_estimation.R
This script estimates confident intervals for alpha coeffitients.

## 07_linear_combination.r
This script calculates summary statistics for shared heredity.

## 07a_validity_check.R
This script provides analytical estimation of heritability of shared heredity and its genetic correlations with original traits using core functions.

## 07b_validity_check.sh
This script estimates heredity of shared heritability and its genetic correlations with original traits using LD Score regression implemented in GWAS-Map database.

## 08_start_gc_correction.sh
This script sources p_correction_for_gc.sh script and sets command line variables (original traits and shared heredity are mentioned).

## 09_p_val_gc_correction.R
This script makes correction for genomic control for original traits and shared heredity.

## 10_clumping.R
This script makes joint clumping for original traits and shared heredity.

## 11_grep_clumped_SNPs.sh
This script extracts summary statistics for clumped SNPs from each of original traits and shared heredity.

## 12_final_clump_table.R
This script forms full table for clumped SNPs with summary statistics from each trait and shared heredity.

## 13_fisher_test_quality_of_SH.R
This script contains Fisher's test estimating quality of shared heredity.

## 14_shared_hits.R
This script looks for shared SNPs for original traits. 

## 15_linear_combination.R
This script substracts shared heredity from original traits.

## 15a_validity_check.R
This script provides analytical estimation of heritability of traits without shared heredity and their genetic correlations with shared heredity using core functions.

## 15b_validity_check.sh
This script estimates heritability of traits without shared heredity and their genetic correlations with shared heredity using LD Score regression implemented in GWAS-Map.

## 16a_start.sh
This script sources 16b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (original traits, shared heredity and traits after substraction of shared heredity are mentioned).

## 16b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits, shared heredity and traits after substraction of shared heredity (see steps 01 - 04). The script uses command line variables.

## 16c_matrix_visualisation.R
This script creates a heatmap plot of genetic correlations of original traits, shared heredity and traits after substraction of shared heredity.

## p_correction_for_gc.sh
This script obtains parameters for genomic control correction for traits set as command line variables. The script is written for GWAS-Map database.

## 17_joint_clumping_and_enrichment.R
This script runs clumping and find shared hits using new algorithm. This script also generates ROC-curve and calculates AUC-value. 
