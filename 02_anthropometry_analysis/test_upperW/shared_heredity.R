source("../../00_core_functions/shared_heredity.R")
path<-'~/polyomica/projects/shared_heredity/data/01_anthropometry_results/four_traits/'
CorPhenTr <- as.matrix(read.table(paste0(path,'pheno_corr_matrix.txt'), check.names=F))
A0 <- as.matrix(read.table(paste0(path,'gene_cov_matrix.txt'), check.names=F))
for (upperW in c(0.6,0.7,0.8,0.9,1)){
	
	output_alpha<-paste0(path, 'upperW_variation/upperW_',upperW,'_alphas.txt')
	output_w<-paste0(path,'upperW_variation/upperW_',upperW,'_weights.txt')
	x<-shared_heredity(CovGenTr=A0, CorPhenTr=CorPhenTr, UpperW=upperW)
	print(x)  
	write.table(x$alphas,output_alpha,quote=F)
	write.table(x$weights,output_w,quote=F)
}
