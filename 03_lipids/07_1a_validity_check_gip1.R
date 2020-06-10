# Aim of this script is to estimate heritability of
# GIP1 and its genetic correlations with other traits

library(data.table)

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/cor_g_a.R")

aa <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/GIP1.txt', row.names = 1) # GIP1 coeffitients
phem <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)
gcov <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

alphas <- as.numeric(aa[, 1])

# Estimate heritability
H2(alphas, covm = gcov, phem = phem)
## 0.2756954

# Estimate pairwise genetic correlations for GIP1 and Lipid traits
cor_gi_alfa(a = alphas, i = 1, covm = gcov) # gip1 and LDL
## 0.9202347
cor_gi_alfa(a = alphas, i = 2, covm = gcov) # gip1 and triglycerides
## 0.736824
cor_gi_alfa(a = alphas, i = 3, covm = gcov) # gip1 and cholesterol
## 0.9154629


