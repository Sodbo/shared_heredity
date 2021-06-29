This folder contains all code for analysis of shared genetic background of lipid traits.
Difference from anthropometric and psychometric traits analysis is in the absence of genomic control p-value correction due to the intercept values <1.

## 00a_start.sh
This script sources 00b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (only original traits are mentioned). 

## 00b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits and estimation of alpha coeffitients for linear combination (see steps 01 - 06). The script uses command line variables.

## 01_pheno_corr.sh
This script calculates phenotypic correlations between traits. The script is writted for GWAS-MAP database and it uses command line variables.

## 02_convert_long_to_wide_form.R
This script converts phenotypic correlations matrix from long to wide format. The script uses command line variables.

## 03_gene_corr.sh
This script calculates genetic correlations between traits and heritability coefficents. Heritability coefficients are taken from ldscroe option --h2 rather than --rg. It is due to the differences in --rg values from both --h2 values and estimated values of heritability. The script is writted for GWAS-MAP database and it uses command line variables.

## 04_gene_corr_to_matrices.R
This script calculates genetic correlations, covariances and se of genetic covariances between traits and builds matrices. The script uses command line variables.

## 05_alpha_coefficients.R
This script calculates alpha coeffitients and weights for each trait in linear combination. The script uses command line variables.

## 06_alpha_CI_estimation.R
This script estimates confident intervals for alpha coeffitients. The script uses command line variables.

## 07_0_linear_combination_sh.R
This script calculates summary statistics for SGCT.

## 07_0a_validity_check_sh.R
This script provides analytical estimation of heritability of SGCT and its genetic correlations with original traits using core functions.

## 07_0b_validity_check_sh.sh
This script estimates heritability of SGCT and its genetic correlations with original traits using LD Score regression implemented in GWAS-MAP database.

## 07_1_linear_combination_gip1.R
This script calculates summary statistics for GIP1 (the GIP approach is described here doi: 10.1038/s42003-020-1051-9).

## 07_1a_validity_check_gip1.R
This script provides analytical estimation of heritability of GIP1 and its genetic correlations with original traits using core functions.

## 07_1b_validity_check_gip1.sh
This script estimates heritability of GIP1 and its genetic correlations with original traits using LD Score regression implemented in GWAS-MAP database.

## 07_1c_rg_sh_gip1.R
This script estimates genetic correlation between GIP1 and SGCT using core function.

## 07_1d_rg_sh_gip1.sh
This script estimates genetic correlation between GIP1 and SGCT using LD Score regression implemented in GWAS-MAP database.

## 07_2_linear_combination_maxh.R
This script calculates summary statistics for MaxH (doi: 10.1159/000381641).

## 07_2a_validity_check_maxh.R
This script provides analytical estimation of heritability of MaxH and its genetic correlations with original traits using core functions.

## 07_2b_validity_check_maxh.sh
This script estimates heritability of MaxH and its genetic correlations with original traits using LD Score regression implemented in GWAS-MAP database.

## 07_2c_rg_sh_maxh.R
This script estimates genetic correlation between MaxH and SGCT using core function.

## 07_2d_rg_sh_maxh.sh
This script estimates genetic correlation between MaxH and SGCT using LD Score regression implemented in GWAS-MAP database.

## 08_linear_combination.R
This script adjusts original traits for SGCT.

## 08a_validity_check.R
This script provides analytical estimation of heritability of UGCTs and their genetic correlations with SGCT using core functions.

## 08b_validity_check.sh
This script estimates heritability of UGCTs and their genetic correlations with SGCT using LD Score regression implemented in GWAS-MAP.

## 09a_start.sh
This script sources 09b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (original traits, SGCT and UGCTs are mentioned).

## 09b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits, SGCT and UGCTs (see steps 01 - 04). The script uses command line variables.

## 09c_matrix_visualisation.R
This script creates a heatmap plot of genetic correlations of original traits, SGCT and UGCTs.

## 10_0_joint_clumping_and_enrichment_sh.R
This script runs clumping for original traits and SGCT and finds shared hits.

## 10_1_joint_clumping_and_enrichment_gip1.R
This script runs clumping for original traits and GIP1 and finds shared hits.

## 10_2_joint_clumping_and_enrichment_maxh.R
This script runs clumping for original traits and MaxH and finds shared hits.

## 11a_count_significant_loci_sh.R
This script counts loci genome-wide significant for each original trait or SGCT.

## 11b_significant_loci_sh.R
This script counts loci genome-wide significant both for particular original trait and SGCT.

## 12_start_gc_correction.sh
This script sources p_correction_for_gc.sh script and sets command line variables (original traits and SGCT are mentioned).

## p_correction_for_gc.sh
This script obtains parameters for genomic control correction for traits set as command line variables. The script is written for GWAS-MAP database. 

## DEPICT
This folder contains config files and script for DEPICT analysis for all lipid traits.
