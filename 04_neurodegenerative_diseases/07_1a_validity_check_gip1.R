# Aim of this script is to estimate heritability of
# GIP1 and its genetic correlations with other traits

library(data.table)

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases")

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/cor_g_a.R")

aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/GIP1.txt', row.names = 1) # GIP1 coeffitients
phem <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)
gcov <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

alphas <- as.numeric(aa[, 1])

# Estimate heritability
H2(alphas, covm = gcov, phem = phem)
## 0.4221662

# Estimate pairwise genetic correlations for GIP1 and PGC traits
cor_gi_alfa(a = alphas, i = 1, covm = gcov) # gip1 and bip
## 0.9486279
cor_gi_alfa(a = alphas, i = 2, covm = gcov) # gip1 and mdd
## 0.4750844
cor_gi_alfa(a = alphas, i = 3, covm = gcov) # gip1 and scz
## 0.8921665
cor_gi_alfa(a = alphas, i = 4, covm = gcov) # gip1 and hap
## 0.3279969

