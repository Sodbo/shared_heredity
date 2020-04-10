# Aim of this script is to estimate correlations between
# shared heredity and the original traits

library(data.table)

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/cor_g_a.R")

aa <- read.table('../../data/01_anthropometry_results/alphas.txt', row.names = 1)

gcov <- read.table('../../data/01_anthropometry_results/gene_cov_matrix.txt', row.names = 1,  check.names = F)

alphas <- as.numeric(aa[2, ])


# Estimate pairwise genetic correlations for SH and anthropometric traits
cor_gi_alfa(a = alphas, i = 1, covm = gcov) # sh and BMI

cor_gi_alfa(a = alphas, i = 2, covm = gcov) # sh and Weight

cor_gi_alfa(a = alphas, i = 3, covm = gcov) # sh and Hip

cor_gi_alfa(a = alphas, i = 4, covm = gcov) # sh and Waist 
