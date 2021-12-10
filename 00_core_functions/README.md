This folder contains core (basic) functions needed to perform analysis using SHAHER framework (including validity checks) and some general functions used specifically for real data analysis.  

The function, that implements MaxSH and GIP aproaches, is shared_heredity.R (see bellow).

The function, that implements sumCOT approch, is linear_combination_v3.R (see bellow).

## cor_g_a.R
This function counts genetic correlation of linear combination with a coefficient "a" with trait "i".

## function_for_estimation_of_alfa_CI.R
This is a function to calculate confidence interval for alpha coeficients.

## gcor_a1_a2.R
This function calculates analytical genetic correlation of linear combination with a coefficient "a1" and a coefficient "a2".

## gcov_a1_a2.R
This function calculates analytical genetic covariance of linear combination with a coefficient "a1" and a coefficient "a2".

## gcov_for_linear_comb_with_i_trait.R
This function counts analytical genetic covariance of linear combination with a coefficient "a" with trait "i".

## heritability_of_linear_combination.R
This function counts analytical heretability of linear combination.

## joint_function_for_enrichment_and_auc.R
This script contains three clumping functions solving different issues:
#### function_for_shlop_29_03_2020
The main function for clumping. Generates a table with clumping results for a set of traits. If several traits from the set contain SNPs, that are located close to each other on a chromosome, only one of this traits with the smallest p-value of association will be presented in the resulting table in the "trait" column. The whole list of such traits will be presented in the "traits" column and the number of these traits will be shown in the "Ntraits" column.
#### clumping_part_I
Part I: clumping of original traits based on given threshold and comparison with SGIT. Calculates the number of loci significantly associated with different number of the original traits and construct a box-plot with -log10(p-value on SGIT) distribution.
#### clumping_part_II
Part II: joint clumping of all traits under p-value < 5e-08 threshold.


## linear_combination_v3.R 
This is a function for sumCOT approach. The function GWAS_linear_combination_Z_based uses Z-scores of the original traits. There is also an alternative function GWAS_linear_combination_v2, which uses beta and se of beta.

## linear_combination_v4.cpp
That is a faster version of GWAS_linear_combination_Z_based function from linear_combination_v3.R script implemented in C++.

## linear_combination_v4_example.R
That is an example of R script which utilizes linear_combination_v4.cpp code for faster sumCOT calculations. 

## p_correction_for_gc.sh
This function estimates SNP-based heritability of the trait and assesses intercept parameter using LD Score regression tool implemented into GWAS-MAP platform (doi: 10.18699/VJ20.686). The function is written for GWAS-MAP platform.

## run_tests.sh
The function for running R scripts from "tests" folder.

## shared_heredity.R
This function estimates weights and tests whether there is a SGI for a set of original traits. If it exists the function estimates alpha coefficients of SGIT linear combination using MaxSH approach. The function also calculates alpha coefficients for MaxH (doi: 10.1159/000381641) and GIP1 (doi: 10.1038/s42003-020-1051-9) traits. The function uses matrices of genetic covariances and phenotypic correlations between the original traits as an input. 

## tests
This folder contains scripts testing the shared_heredity.R function for errors and reproducibility of the results.
