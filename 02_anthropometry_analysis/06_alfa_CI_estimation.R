# Aim of this script is estimation of confidence intervals for alpha coefficients in linear combination

#gcovm - matrix of genetics covariances
#phem - phneotypic covs
#se - matrix of se of h2 and gcovm
#lll - shrinkage factor [0;1]
#N_permut - N of permutations

path<-'../../data/01_anthropometry_results/five_traits/'
CorPhenTr <- as.matrix(read.table(paste0(path,'pheno_corr_matrix.txt'), check.names=F))
A0 <- as.matrix(read.table(paste0(path,'gene_cov_matrix.txt'), check.names=F))

source('../00_core_functions/shared_heredity.R', chdir = F)
se <-as.matrix(read.table(paste0(path,'gene_cov_se_matrix.txt'), check.names=F))

source('../00_core_functions/function_for_estimation_of_alfa_CI.R', chdir = F)
  (res<-function_for_estimation_of_alfa_CI(A0,CorPhenTr,se,N_permut = 100))
write.table(res,paste0(path,'CIs_for_5_traits.txt'),quote=F)

