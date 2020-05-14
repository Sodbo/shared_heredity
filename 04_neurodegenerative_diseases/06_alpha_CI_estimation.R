# Aim of this script is estimation of confidence intervals for alpha coefficients in linear combination

#gcovm - matrix of genetics covariances
#phem - phneotypic covs
#se - matrix of se of h2 and gcovm
#lll - shrinkage factor [0;1]
#N_permut - N of permutations

path <- '/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/four_traits/'

CorPhenTr <- as.matrix(read.table(paste0(path, 'pheno_corr_matrix.txt'), check.names = F)) # load matrix of phenotypic correlations
A0 <- as.matrix(read.table(paste0(path, 'gene_cov_matrix.txt'), check.names = F)) # load matrix of genetic covariance

se <- as.matrix(read.table(paste0(path,'gene_cov_se_matrix.txt', check.names = F))

source("../00_core_functions/function_for_estimation_of_alfa_CI.R")

res <- function_for_estimation_of_alfa_CI(gcovm = A0, phem = CorPhenTr, se = se, N_permut = 100)
names(res) <- colnames(A0)
write.table(res, '/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/four_traits/CIs_for_4_traits.txt', quote = F)

