# Aim of this script is to estimate heritability of traits after
# shared heredity subtraction and their genetic correlations with shared heredity

library(data.table)


source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/gcor_a1_a2.R")
source("../00_core_functions/gcov_for_linear_comb_with_i_trait.R")

aa <- read.table('../../data/01_anthropometry_results/alphas.txt', row.names = 1)

alphas <- as.numeric(aa[2, ])

n_traits <- c(1:length(alphas))

phem <- read.table('../../data/01_anthropometry_results/pheno_corr_matrix.txt', row.names = 1,  check.names = F)

gcov <- read.table('../../data/01_anthropometry_results/gene_cov_matrix.txt', row.names = 1,  check.names = F)

slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))

position <- diag(length(alphas))

# Estimate heritability
tmp <- lapply(n_traits, function(x) H2(position[x, ] - alphas*slope[x], covm = gcov, phem = phem))
as.numeric(tmp)
#0.15165770 0.12550849 0.08064241 0.07232553
# Estimate pairwise genetic correlations for traits-SH and SH
tmp2 <- lapply(n_traits, function(x) cor_gi_a1_a2(a1 = alphas, a2 = position[x, ] - alphas*slope[x], covm = gcov))
as.numeric(tmp2)
