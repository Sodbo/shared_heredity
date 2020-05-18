#gcovm - matrix of genetics covariances
#phem - phneotypic covs
#se - matrix of se of h2 and gcovm
#lll - shrinkage factor [0;1]
#N_permut - N of permutations

source("05_shared_heredity.R", chdir = F)
se <-as.matrix(read.table('~/polyomica/projects/shared_heredity/data/02_Lipids/gene_cov_se_matrix.txt', check.names=F))
#gcovm <-read.table('~/polyomica/projects/shared_heredity/data/02_Lipids/gene_cov_matrix.txt', check.names=F)
#phem <-read.table('~/polyomica/projects/shared_heredity/data/02_Lipids/pheno_corr_matrix.txt', check.names=F)

source("../00_core_functions/function_for_estimation_of_alfa_CI.R")

res<-function_for_estimation_of_alfa_CI(A0,CorPhenTr,se,N_permut = 200)
names(res)=colnames(A0)
write.table(res,'~/polyomica/projects/shared_heredity/data/02_Lipids/CIs_for_3_traits.txt',quote=F)

