# Aim of this script is to estimate heritability of
# MaxH and its genetic correlations with other traits

library(data.table)

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/cor_g_a.R")

aa <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/alphas.txt', row.names = 1)
phem <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)
gcov <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

alphas <- as.numeric(aa[1, ])

# Estimate heritability
H2(alphas, covm = gcov, phem = phem)
## 0.3011047

# Estimate pairwise genetic correlations for MaxH and PGC traits
cor_gi_alfa(a = alphas, i = 1, covm = gcov) # maxh and LDL
## 0.7983689
cor_gi_alfa(a = alphas, i = 2, covm = gcov) # maxh and triglycerides
## 0.8887436
cor_gi_alfa(a = alphas, i = 3, covm = gcov) # maxh and cholesterol
## 0.7846558

