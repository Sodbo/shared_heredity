#!/bin/bash
# Path to output directory should be given as the first command line argument
# Ids of GWAS should also be given as the command line arguments after path.
# Necessary R libraries should be installed before the first run.
# Results will be written to the following files:
# gene_corr_matrix.txt
# gene_corr directory
# gene_cov_matrix.txt
# phen_corr_res.txt
# pheno_corr_matrix.txt
# alphas.txt
# w.txt

## Example
## source 14b_pipeline_for_matrix_calculation.sh ../../data/anthropometry_results/four_traits/Traits_minus_SH_test/ 153 154 155 156

source pipeline_for_calculation_of_matrices.sh ../../data/01_anthropometry_results/five_traits/Traits_vs_Traits_minus_SH/ 191 192 193 194 199 201 202 203 204 205 206 

