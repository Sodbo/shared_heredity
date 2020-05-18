# Aim of this script is to estimate heritability of traits after
# shared heredity subtraction and their genetic correlations with shared heredity

library(data.table)

source("../00_core_functions/heritability_of_linear_combination.R")
source("../00_core_functions/gcor_a1_a2.R")
source("../00_core_functions/gcov_for_linear_comb_with_i_trait.R")

aa <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/alphas.txt', row.names = 1)

alphas <- as.numeric(aa[2, ])

n_traits <- c(1:length(alphas))

phem <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)

gcov <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))

position <- diag(length(alphas))

# Estimate heritability
tmp <- lapply(n_traits, function(x) H2(position[x, ] - alphas*slope[x], covm = gcov, phem = phem))
#[[1]]
#[1] 0.0863137 for LDL-sh
#[[2]]
#[1] 0.1945357 for triglycerides-sh
#[[3]]
#[1] 0.1698012 for cholesterol-sh

# Estimate pairwise genetic correlations for traits-SH and SH
tmp2 <- lapply(n_traits, function(x) cor_gi_a1_a2(a1 = alphas, a2 = position[x, ] - alphas*slope[x], covm = gcov))
# [[1]]
# [1] 1.741672e-14 rg for sh and LDL-sh
# [[2]]
# [1] -2.819108e-15 rg for sh and triglycerides-sh
# [[3]]
# [1] -1.349233e-14 rg for sh and cholesterol-sh

