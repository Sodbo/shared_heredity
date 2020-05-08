source("../00_core_functions/shared_heredity.R")
path<-commandArgs(trailingOnly=T)[1]
CorPhenTr <- as.matrix(read.table(paste0(path,'pheno_corr_matrix.txt'), check.names=F))
A0 <- as.matrix(read.table(paste0(path,'gene_cov_matrix.txt'), check.names=F))
output_alpha<-paste0(path,'alphas.txt')
output_w<-paste0(path,'weights.txt')
output_gip<-paste0(path, 'GIP1.txt')


 x<-shared_heredity(CovGenTr=A0, CorPhenTr=CorPhenTr)
print(x)  
write.table(x$alphas,output_alpha,quote=F)
write.table(x$weights,output_w,quote=F)
write.table(x$GIPs$GIP_coeff[,'GIP1'], output_gip, quote=F)
x$alphas
