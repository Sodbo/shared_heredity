source("../00_core_functions/shared_heredity.R")
CorPhenTr <- as.matrix(read.table('~/polyomica/projects/shared_heredity/data/02_Lipids/pheno_corr_matrix.txt', check.names=F))
A0 <- as.matrix(read.table('~/polyomica/projects/shared_heredity/data/02_Lipids/gene_cov_matrix.txt', check.names=F))

x<-shared_heredity(CovGenTr=A0, CorPhenTr=CorPhenTr)
  
write.table(x$alphas,'~/polyomica/projects/shared_heredity/data/02_Lipids/alphas.txt',quote=F)
write.table(x$weights,'~/polyomica/projects/shared_heredity/data/02_Lipids/weights.txt',quote=F)

