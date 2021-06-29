# Aim of this script is to estimate heritability of
# MaxH and its genetic correlations with other traits

library(data.table)

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases")

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/cor_g_a.R")

aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/alphas.txt', row.names = 1)
phem <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)
gcov <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

alphas <- as.numeric(aa[1, ])

# Estimate heritability
H2(alphas, covm = gcov, phem = phem)
## 0.4256747

# Estimate pairwise genetic correlations for MaxH and PGC traits
cor_gi_alfa(a = alphas, i = 1, covm = gcov) # maxh and bip
## 0.9566709
cor_gi_alfa(a = alphas, i = 2, covm = gcov) # maxh and mdd
## 0.4565465
cor_gi_alfa(a = alphas, i = 3, covm = gcov) # maxh and scz
## 0.8831441
cor_gi_alfa(a = alphas, i = 4, covm = gcov) # maxh and hap
## 0.3128241

