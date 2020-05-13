# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with other traits

library(data.table)

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases")

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/cor_g_a.R")

aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/alphas.txt', row.names = 1)
phem <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)
gcov <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

alphas <- as.numeric(aa[2, ])

# Estimate heritability
H2(alphas, covm = gcov, phem = phem)
## 0.4167322 (4 traits sh)

# Estimate pairwise genetic correlations for SH and PGC traits
cor_gi_alfa(a = alphas, i = 1, covm = gcov) # sh and bip
## 0.9414797 (4 traits sh)
cor_gi_alfa(a = alphas, i = 2, covm = gcov) # sh and mdd
## 0.495251 (4 traits sh)
cor_gi_alfa(a = alphas, i = 3, covm = gcov) # sh and scz
## 0.8956868 (4 traits sh)
cor_gi_alfa(a = alphas, i = 4, covm = gcov) # sh and hap
## 0.3573614

