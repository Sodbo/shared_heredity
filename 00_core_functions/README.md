This folder contains core (basic) functions needed to perform analysis using MaSH framework (including validity checks) and some general functions used specifically for real data analysis.  

## cor_g_a.R
This function counts genetic correlation of linear combination with a coefficient "a" with trait "i".

## function_for_estimation_of_alfa_CI.R
This is a function to calculate confidence interval for alpha coeficients.

## gcor_a1_a2.R
This function calculates analytical genetic correlation of linear combination with a coefficient "a1" and a coefficient "a2".

## gcov_a1_a2.R
This function calculates analytical genetic covariation of linear combination with a coefficient "a1" and a coefficient "a2".

## gcov_for_linear_comb_with_i_trait.R
This function counts analytical genetic covariance of linear combination with a coefficient "a" with trait "i".

## heritability_of_linear_combination.R
This function counts analytical heretability of linear combination.

## joint_function_for_enrichment_and_auc.R
This script contains three clumping functions solving different issues:
### function_for_shlop_29_03_2020
The main function for clumping. Generates a table with clumping results for a set of traits. If several traits from the set contain SNPs, that are located close to each other on a chromosome, only one of this traits with the smallest p-value of association will be presented in the resulting table in the "trait" column. The whole list of such traits will be presented in the "traits" column and the number of these traits will be shown in the "Ntraits" column.
### clumping_part_I
Part I: clumping of original traits based on given threshold and comparison with SGCT. Calculates number of shared hits by N_at level, and construct box-plot with shared hits distribution.
### clumping_part_II
Part II: joint clumping of all traits under p-value < 5e-08 threshold.

## linear_combination_v3.R 
This is a function calculating GWAS summary statistics for SGCT of the original traits based on their summary statistics. The function uses Z-scores of the original traits.

## p_correction_for_gc.sh
This function estimates SNP-based heritability of the trait and assesses intercept parameter using LD Score regression tool implemented into GWAS-MAP platform (doi: 10.18699/VJ20.686). The function is written for GWAS-MAP platform.

## run_tests.sh
The function for running R scripts from "tests" folder.

## shared_heredity.R
This function estimates weights and tests whether there is a SGC for a set of original traits. If it exists the function estimates alpha coefficients of SGCT linear combination. The function also calculates alpha coefficients for MaxH (doi: 10.1159/000381641) and GIP1 (doi: 10.1038/s42003-020-1051-9) traits. The function uses matrices of genetic covariances and phenotypic correlations between the original traits as an input. 

## tests
This folder contains scripts testing the shared_heredity.R function for errors and reproducibility of the results.
