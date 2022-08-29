This folder contains all code for analysis of shared genetic background of psychometric traits. Note: SH is an abbreviation of "shared heredity". Used in the code as the synonym of SGIT.

## 00a_start.sh
This script sources 00b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (only original traits are mentioned). 

## 00b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits and estimation of alpha coeffitients for linear combination (see steps 01 - 05). The script uses command line variables.

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

## 06_alpha_CI_estimation.R
This script estimates confident intervals for alpha coeffitients.

## 07_0_linear_combination_sh.R
This script calculates summary statistics for SGIT.

## 07_0a_validity_check_sh.R
This script provides analytical estimation of heritability of SGIT and its genetic correlations with original traits using core functions.

## 07_0b_validity_check_sh.sh
This script estimates heritability of SGIT and its genetic correlations with original traits using LD Score regression implemented in GWAS-MAP database.

## 07_1_linear_combination_gip1.R
This script calculates summary statistics for GIP1 (the GIP approach is described here doi: 10.1038/s42003-020-1051-9).

## 07_1a_validity_check_gip1.R
This script provides analytical estimation of heritability of GIP1 and its genetic correlations with original traits using core functions.

## 07_1b_validity_check_gip1.sh
This script estimates heritability of GIP1 and its genetic correlations with original traits using LD Score regression implemented in GWAS-MAP database.

## 07_1c_rg_sh_gip1.R
This script estimates genetic correlation between GIP1 and SGIT using core function.

## 07_1d_rg_sh_gip1.sh
This script estimates genetic correlation between GIP1 and SGIT using LD Score regression implemented in GWAS-MAP database.

## 07_2_linear_combination_maxh.R
This script calculates summary statistics for MaxH (doi: 10.1159/000381641).

## 07_2a_validity_check_maxh.R
This script provides analytical estimation of heritability of MaxH and its genetic correlations with original traits using core functions.

## 07_2b_validity_check_maxh.sh
This script estimates heritability of MaxH and its genetic correlations with original traits using LD Score regression implemented in GWAS-MAP database.

## 07_2c_rg_sh_maxh.R
This script estimates genetic correlation between MaxH and SGIT using core function.

## 07_2d_rg_sh_maxh.sh
This script estimates genetic correlation between MaxH and SGIT using LD Score regression implemented in GWAS-MAP database.

## 08_0_start_gc_correction.sh
This script sources p_correction_for_gc.sh script and sets command line variables (original traits and SGIT are mentioned).

## 08_1_start_gc_correction_gip1.sh
This script sources p_correction_for_gc.sh script and sets command line variables (GIP1 is mentioned).

## 08_2_start_gc_correction_maxh.sh
This script sources p_correction_for_gc.sh script and sets command line variables (MaxH is mentioned).

## 09_0_p_val_gc_correction.R
This script makes correction for genomic control for original traits, SGIT and UGITs.

## 09_1_p_val_gc_correction_gip1.R
This script makes correction for genomic control for GIP1.

## 09_2_p_val_gc_correction_maxh.R
This script makes correction for genomic control for MaxH.

## 10_linear_combination.R
This script adjusts original traits for SGIT.

## 10a_validity_check.R
This script provides analytical estimation of heritability of traits adjusted for SGIT and their genetic correlations with SGIT using core functions.

## 10b_validity_check.sh
This script estimates heritability of traits adjusted for SGIT and their genetic correlations with SGIT using LD Score regression implemented in GWAS-MAP.

## 11a_start.sh
This script sources 11b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (original traits, SGIT and UGITs are mentioned).

## 11b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits, SGIT and UGITs (see steps 01 - 04). The script uses command line variables.

## 11c_matrix_visualisation.R
This script creates a heatmap plot of genetic correlations of original traits, SGIT and UGITs.

## p_correction_for_gc.sh
This script obtains parameters for genomic control correction for traits set as command line variables. The script is written for GWAS-MAP database.

## 12_0_joint_clumping_and_enrichment_sh.R
This script runs clumping for original traits and SGIT and finds loci significantly associated fith different number of original traits.

## 12a_0_joint_clumping_and_enrichment_sh.R
This script runs clumping for original traits, SGIT and UGITs and finds loci significantly associated fith different number of original traits.

## 12_1_joint_clumping_and_enrichment_gip1.R
This script runs clumping for original traits and GIP1 and finds loci significantly associated fith different number of original traits.

## 12_2_joint_clumping_and_enrichment_maxh.R
This script runs clumping for original traits and MaxH and finds loci significantly associated fith different number of original traits.

## 13_count_significant_loci.R
This script counts genome-wide significant loci for each original trait or SGIT.

## 14_significant_loci_sh.R
This script counts loci genome-wide significant both for particular original trait and for SGIT.

## 15_matrix_visualisation_shared_loci.R
This script performs a heatmap visualization of loci overlap between traits. 

## 16_check_new_snps_in_bigger_gwas.R
This script extracts p-values of SNPs genowide-significantly associated with at least one of the traits from the psychometry dataset from the bigger GWAS on MDD.

## DEPICT
This folder contains config files and script for DEPICT analysis for all psychometric traits.
