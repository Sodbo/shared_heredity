# Aim of this script is estimation of confidence intervals for alpha coefficients in linear combination

#gcovm - matrix of genetics covariances
#phem - phneotypic covs
#se - matrix of se of h2 and gcovm
#lll - shrinkage factor [0;1]
#N_permut - N of permutations

path <- commandArgs(trailingOnly = T)[1] # get command arguments set in 00a_start.sh or 16a_start.sh

CorPhenTr <- as.matrix(read.table(paste0(path, 'pheno_corr_matrix.txt'), check.names = F)) # load matrix of phenotypic correlations
A0 <- as.matrix(read.table(paste0(path, 'gene_cov_matrix.txt'), check.names = F)) # load matrix of genetic covariance

se <- as.matrix(read.table(paste0(path,'gene_cov_se_matrix.txt'), check.names = F))

source("../00_core_functions/function_for_estimation_of_alfa_CI.R")
source("../00_core_functions/shared_heredity.R")

# for alpha coeffitients of SH
res <- function_for_estimation_of_alfa_CI(gcovm = A0, phem = CorPhenTr, se = se, N_permut = 100)
names(res) <- colnames(A0)
write.table(res, paste0(path,'CIs_for_',length(A0[1,]),'_traits.txt'), quote = F)

# for alpha coeffitients of GIP1
res2 <- function_for_estimation_of_alfa_CI_GIP(gcovm = A0, phem = CorPhenTr, se = se, N_permut = 100)
names(res2) <- colnames(A0)
ifelse(!dir.exists(file.path(path, 'GIP')), dir.create(file.path(path, 'GIP')), FALSE)
write.table(res2, paste0(path, 'GIP/CIs_for_',length(A0[1,]),'_traits.txt'), quote = F)

# for alpha coeffitients of MaxH
res3 <- function_for_estimation_of_alfa_CI_MaxH(gcovm = A0, phem = CorPhenTr, se = se, N_permut = 100)
names(res3) <- colnames(A0)
ifelse(!dir.exists(file.path(path, 'MaxH')), dir.create(file.path(path, 'MaxH')), FALSE)
write.table(res3, paste0(path,'MaxH/CIs_for_',length(A0[1,]),'_traits.txt'), quote = F)

