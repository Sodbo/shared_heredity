
### cor_g_a.R
Function counts genetic correlation of linear combination with a coefficient "a" with trait "i"

###function_for_estimation_of_alfa_CI.R
Function to calculate confidence interval for alpha coeficients

###gcor_a1_a2.R
Function calculates genetic correlation of linear combination with a coefficient "a1" and a coefficient "a2"

###gcov_a1_a2.R
Function calculates genetic covariation of linear combination with a coefficient "a1" and a coefficient "a2"

###gcov_for_linear_comb_with_i_trait.R
Counts  genetic covariance of linear combination with a coefficient "a" with trait "i"

###heritability_of_linear_combination.R
Counts heretability of linear combination

###joint_function_for_enrichment_and_auc.R
Contains three function of clumping for different tasks:
## function_for_shlop_29_03_2020
Main function for clumping. If you apply this function with parameter "trait" it takes its value as a column name with trait designation to clump throughout several traits. If two traits contained SNPs, that are close to each other on chromosome, only one of them will be presented in result table and the smallest p-value. But both traits will be shown in "traits" column and the best one will be shown in "trait" column. Also "Ntraits" will contain the number of traits with the close SNPs.
## clumping_part_I
Part I: clumping of original traits based on given threshold and comparison with SGCT. Calculates number of shared hits by N_at level, and construct box-plot with shared hit distribution.

## clumping_part_II
Part II: joint clumping of all traits on 5e-8 threshold

### linear_combination_v3.R 
A function to construct a GWAS summary statistics file for the general genetic component of analyzed phenotypic traits as a linear combination of beta values from summary statistics files for these individual phenotypic traits. This version using as input parameter Z-scores of the analysed traits.

###p_correction_for_gc.sh


###run_tests.sh

###shared_heredity.R

###tests
