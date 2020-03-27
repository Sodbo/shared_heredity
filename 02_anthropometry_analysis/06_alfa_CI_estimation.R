#gcovm - matrix of genetics covariances
#phem - phneotypic covs
#se - matrix of se of h2 and gcovm
#lll - shrinkage factor [0;1]
#N_permut - N of permutations

source("05_shared_heredity.R", chdir = F)
se <-as.matrix(read.table('../data/anthropometry_results/gene_cov_se_matrix.txt', check.names=F))

source('../00_core_functions/gcov_for_linear_comb_with_i_trait.R', chdir = F)
  (res<-function_for_estimation_of_alfa_CI(A0,CorPhenTr,se,N_permut = 3))
    wanames(res)<-colnames()
write.table(res,'../data/anthropometry_results/CIs_for_5_traits.txt',quote=F)

