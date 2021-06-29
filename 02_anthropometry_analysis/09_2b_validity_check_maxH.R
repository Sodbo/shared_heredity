# Aim of this script is to estimate correlations between
# maximal heritability component and the original traits

library(data.table)

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/cor_g_a.R")

aa <- read.table('../../data/01_anthropometry_results/five_traits/alphas.txt', row.names = 1)

gcov <- read.table('../../data/01_anthropometry_results/five_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

alphas <- as.numeric(aa[1, ])


# Estimate pairwise genetic correlations for maxH and anthropometric traits
cor_gi_alfa(a = alphas, i = 1, covm = gcov) # maxh and BMI

cor_gi_alfa(a = alphas, i = 2, covm = gcov) # maxh and Weight

cor_gi_alfa(a = alphas, i = 3, covm = gcov) # maxh and Hip

cor_gi_alfa(a = alphas, i = 4, covm = gcov) # maxh and Waist 

cor_gi_alfa(a = alphas, i = 5, covm = gcov) # maxh and Fat 
