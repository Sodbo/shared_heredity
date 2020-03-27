# Aim of this script is estimation of confidence intervals for alpha coefficients in linear combination

#gcovm - matrix of genetics covariances
#phem - phneotypic covs
#se - matrix of se of h2 and gcovm
#lll - shrinkage factor [0;1]
#N_permut - N of permutations

source("05_alpha_coefficients.R", chdir = F) # variables A0 and CorPhenTr are obtained from this script

se <-as.matrix(read.table('/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/gene_cov_se_matrix.txt', check.names = F))

source("../00_core_functions/function_for_estimation_of_alfa_CI.R")

res <- function_for_estimation_of_alfa_CI(gcovm = A0, phem = CorPhenTr, se = se, N_permut = 100)
names(res) <- colnames(A0)
write.table(res, '/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/CIs_for_3_traits.txt', quote = F)

