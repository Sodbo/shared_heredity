### clumping.R
Function for clumping. This is version 2 with simplified code in comparison to old one from commit 10300ae175d92e61a5e66ebd19e8cb9f11832829. It was tested and show the same results as previous one. If you apply this function with parameter "trait" it takes its value as a column name with trait designation to clump throughout several traits. If two traits contained SNPs, that are close to each other on chromosome, only one of them will be presented in result table and the smallest p-value. But both traits will be shown in "traits" column and the best one will be shown in "trait" column. Also "Ntraits" will contain the number of traits with the close SNPs.

### linear_combination.R 
A function to construct a GWAS summary statistics file for the general genetic component of analyzed phenotypic traits as a linear combination of beta values from summary statistics files for these individual phenotypic traits.This version using as input parameters phenotypic variance, effect sizes and standard errors of the analysed traits.

### linear_combination_v3.R 
A function to construct a GWAS summary statistics file for the general genetic component of analyzed phenotypic traits as a linear combination of beta values from summary statistics files for these individual phenotypic traits. This version using as input parameter Z-scores of the analysed traits.


