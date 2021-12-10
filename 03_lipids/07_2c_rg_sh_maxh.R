# Aim of this script is to estimate genetic correlation between
# SGIT and MaxH

library(data.table)

source("../00_core_functions/gcor_a1_a2.R")

aa <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/alphas.txt', row.names = 1)
gcov <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

maxh_coef <- as.numeric(aa[1, ])
alphas <- as.numeric(aa[2, ])

# Estimate pairwise genetic correlation for MaxH and SH
cor_gi_a1_a2(a1 = maxh_coef, a2 = alphas, covm = gcov)
## 0.9269505


