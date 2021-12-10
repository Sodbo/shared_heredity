# Aim of this script is to estimate genetic correlation between
# SGIT and MaxH

library(data.table)

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases")

source("../00_core_functions/gcor_a1_a2.R")

aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/alphas.txt', row.names = 1)
gcov <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

maxh_coef <- as.numeric(aa[1, ])
alphas <- as.numeric(aa[2, ])

# Estimate pairwise genetic correlation for MaxH and SGIT
cor_gi_a1_a2(a1 = maxh_coef, a2 = alphas, covm = gcov)
## 0.9980714


