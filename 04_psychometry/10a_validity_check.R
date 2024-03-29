# Aim of this script is to estimate heritability of UGITs
# and their genetic correlations with SGIT

library(data.table)

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases")

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/gcor_a1_a2.R")
source("../00_core_functions/gcov_for_linear_comb_with_i_trait.R")

aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/alphas.txt', row.names = 1)

alphas <- as.numeric(aa[2, ])

n_traits <- c(1:length(alphas))

phem <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)

gcov <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))

position <- diag(length(alphas))

# Estimate heritability
tmp <- lapply(n_traits, function(x) H2(position[x, ] - alphas*slope[x], covm = gcov, phem = phem))
#[[1]]
#[1] 0.1078178 for bip UGIT
#[[2]]
#[1] 0.05786701 for mdd UGIT
#[[3]]
#[1] 0.0987432 for scz UGIT
#[[4]]
#[1] 0.05639504 for happiness UGIT


# Estimate pairwise genetic correlations for UGITs and SGIT
tmp2 <- lapply(n_traits, function(x) cor_gi_a1_a2(a1 = alphas, a2 = position[x, ] - alphas*slope[x], covm = gcov))
# [[1]]
# [1] 5.312335e-15 rg for SGIT and bip UGIT
# [[2]]
# [1] -8.813521e-16 rg for SGIT and mdd UGIT
# [[3]]
# [1] -4.858323e-15 rg for SGIT and scz UGIT
# [[4]]
# [1] 2.81601e-16 rg for SGIT and happiness UGIT

