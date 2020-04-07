# Aim of this script is to estimate heritability of traits after
# shared heredity subtraction and their genetic correlations with shared heredity

library(data.table)

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases")

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/gcor_a1_a2.R")
source("../00_core_functions/gcov_for_linear_comb_with_i_trait.R")

aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/alphas.txt', row.names = 1)

alphas <- as.numeric(aa[2, ])

n_traits <- c(1:length(alphas))

phem <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)

gcov <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))

position <- diag(length(alphas))

# Estimate heritability
tmp <- lapply(n_traits, function(x) H2(position[x, ] - alphas*slope[x], covm = gcov, phem = phem))
# [[1]]       
# [1] 0.1039611  h2 for bip-sh
# [[2]]
# [1] 0.05936681 h2 for mdd-sh
# [[3]]
# [1] 0.1029796 h2 for scz-sg

# Estimate pairwise genetic correlations for traits-SH and SH
tmp2 <- lapply(n_traits, function(x) cor_gi_a1_a2(a1 = alphas, a2 = position[x, ] - alphas*slope[x], covm = gcov))
# [[1]]
# [1] 5.05388e-15 rg for sh and bip-sh
# [[2]]
# [1] -2.517877e-15 rg for sh and mdd-sh
# [[3]]
# [1] -4.455682e-15 rg for sh and scz-sh

