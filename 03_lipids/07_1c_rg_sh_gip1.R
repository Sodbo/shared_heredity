# Aim of this script is to estimate genetic correlation between
# shared heredity and GIP1

library(data.table)

source("../00_core_functions/gcor_a1_a2.R")

gip1 <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/GIP1.txt', row.names = 1) # GIP1 coeffitients
aa <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/alphas.txt', row.names = 1) # SH coeffitients
gcov <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

gip1_coef <- as.numeric(gip1[, 1])
alphas <- as.numeric(aa[2, ])

# Estimate pairwise genetic correlation for GIP1 and SH
cor_gi_a1_a2(a1 = gip1_coef, a2 = alphas, covm = gcov)
## 0.993006


