This folder contains all code extracting shared heredity from Lipid traits, excepting uploading gwases to container which are in the ../01_upload_GWAS/Lipids directory.
Difference from anthropometry and PGC pipelines is in absence of genomic control p-value correction due to the intercept values <1.

## 00a_start.sh
This script sources 00b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (only original traits are mentioned). 

## 00b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits and estimation of alpha coeffitients for linear combination (see steps 01 - 06). The script uses command line variables.

## 01_pheno_corr.sh
This script calculates phenotypic correlations between traits. The script is writted for GWAS-Map database and it uses command line variables.

## 02_convert_long_to_wide_form.R
This script converts phenotypic correlations matrix from long to wide format. The script uses command line variables.

## 03_gene_corr.sh
This script calculates genetic correlations between traits and heritability coefficents. Heritability coefficients are taken from ldscroe option --h2 rather than --rg. It is due to the differences in --rg values from both --h2 values and estimated values of heritability. The script is writted for GWAS-Map database and it uses command line variables.

## 04_gene_corr_to_matrices.R
This script calculates genetic correlations, covariances and se of genetic covariances between traits and builds matrices. The script uses command line variables.

## 05_alpha_coefficients.R
This script calculates alpha coeffitients and weights for each trait in linear combination. The script uses command line variables.

## 06_alpha_CI_estimation.R
This script estimates confident intervals for alpha coeffitients. The script uses command line variables.

## 07_0_linear_combination_sh.R
This script calculates summary statistics for shared heredity.

## 07_0a_validity_check_sh.R
This script provides analytical estimation of heritability of shared heredity and its genetic correlations with original traits using core functions.

## 07_0b_validity_check_sh.sh
This script estimates heritability of shared heredity and its genetic correlations with original traits using LD Score regression implemented in GWAS-Map database.

## 07_1_linear_combination_gip1.R
This script calculates summary statistics for GIP1.

## 07_1a_validity_check_gip1.R
This script provides analytical estimation of heritability of GIP1 and its genetic correlations with original traits using core functions.

## 07_1b_validity_check_gip1.sh
This script estimates heritability of GIP1 and its genetic correlations with original traits using LD Score regression implemented in GWAS-Map database.

## 07_1c_rg_sh_gip1.R
This script estimates genetic correlation between GIP1 and shared heredity using core function.

## 07_1d_rg_sh_gip1.sh
This script estimates genetic correlation between GIP1 and shared heredity using LD Score regression implemented in GWAS-Map database.

## 07_2_linear_combination_maxh.R
This script calculates summary statistics for MaxH.

## 07_2a_validity_check_maxh.R
This script provides analytical estimation of heritability of MaxH and its genetic correlations with original traits using core functions.

## 07_2b_validity_check_maxh.sh
This script estimates heritability of MaxH and its genetic correlations with original traits using LD Score regression implemented in GWAS-Map database.

## 07_2c_rg_sh_maxh.R
This script estimates genetic correlation between MaxH and shared heredity using core function.

## 07_2d_rg_sh_maxh.sh
This script estimates genetic correlation between MaxH and shared heredity using LD Score regression implemented in GWAS-Map database.

## 08_linear_combination.R
This script substracts shared heredity from original traits.

## 08a_validity_check.R
This script provides analytical estimation of heritability of traits without shared heredity and their genetic correlations with shared heredity using core functions.

## 08b_validity_check.sh
This script estimates heritability of traits without shared heredity and their genetic correlations with shared heredity using LD Score regression implemented in GWAS-Map.

## 09a_start.sh
This script sources 09b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (original traits, shared heredity and traits after substraction of shared heredity are mentioned).

## 09b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits, shared heredity and traits after substraction of shared heredity (see steps 01 - 04). The script uses command line variables.

## 09c_matrix_visualisation.R
This script creates a heatmap plot of genetic correlations of original traits, shared heredity and traits after substraction of shared heredity.

## 10_0_joint_clumping_and_enrichment_sh.R
This script runs clumping for original traits and shared heredity and find shared hits. This script also generates ROC-curve and calculates AUC-value.

## 10_1_joint_clumping_and_enrichment_gip1.R
This script runs clumping for original traits and GIP1 and find shared hits. This script also generates ROC-curve and calculates AUC-value.

## 10_2_joint_clumping_and_enrichment_maxh.R
This script runs clumping for original traits and MaxH and find shared hits. This script also generates ROC-curve and calculates AUC-value.
 
