# Aim of this script is to estimate genetic correlation between
# SGIT and GIP1

library(data.table)

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases")

source("../00_core_functions/gcor_a1_a2.R")

gip1 <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/GIP1.txt', row.names = 1) # GIP1 coeffitients
aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/alphas.txt', row.names = 1) # SGIT coeffitients
gcov <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

gip1_coef <- as.numeric(gip1[, 1])
alphas <- as.numeric(aa[2, ])

# Estimate pairwise genetic correlation for GIP1 and SGIT
cor_gi_a1_a2(a1 = gip1_coef, a2 = alphas, covm = gcov)
## 0.9994053 


