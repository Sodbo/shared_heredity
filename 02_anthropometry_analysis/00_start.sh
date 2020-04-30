#!/bin/bash
# This script starts a pipline script to count matrices of genetic and phenotypic covariations required for shared genetic component and also start the script to calculate it. Check the starting parameters of this script enumerating 05 if necessary.
# Path to output directory should be given as the first command line argument
# Ids of GWAS should also be given as the command line arguments after path.
# Necessary R libraries should be installed before the first run.
# 
# Results will be written to the following files (path to them you should give as the first parameter):
# gene_corr_matrix.txt
# gene_corr directory
# gene_cov_matrix.txt
# phen_corr_res.txt
# pheno_corr_matrix.txt
# alphas.txt
# w.txt

## Example
## source pipeline_for_matrix_calculation.sh ../../data/01_anthropometry_results/Traits_minus_SH_test/ 153 154 155 156
#export PROD=T #Option is necessary for GWASes uploaded to GWAS-Map
source pipeline_for_calculation_of_matrices.sh ../../data/01_anthropometry_results/five_traits/ 191 192 193 194 199

