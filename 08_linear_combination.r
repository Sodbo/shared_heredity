library(data.table)
library(dplyr)

path_to_result_directory<-'../data/anthropometry_results/four_traits/GWAS/scaled_filtered/'
gwas_files<-list.files(path_to_result_directory, full.names=T, pattern="ID_\\d+.csv")
gwas<-lapply(gwas_files, fread)

aa<-read.table('../data/anthropometry_results/four_traits/alphas.txt', row.names=1)
covm<-read.table('../data/anthropometry_results/four_traits/pheno_corr_matrix.txt', row.names=1,  check.names=F)

rs_id<-lapply(gwas, function(x) x$rs_id)
snps<-Reduce(intersect,rs_id)
ind<-lapply(rs_id, function(x) match(snps,x))

#Is it necessary to remove NA for the case of not whole intersection?
#ind<-lapply(ind,function(x) )

gwas_reordered<-lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]],])
betas<-sapply(gwas_reordered, function(x) x$beta)
se <- sapply(gwas_reordered, function(x) x$se)
sample_size <- sapply(gwas_reordered, function(x) x$n)

GWAS_linear_combination=function(a,beta_a,se,var_y=rep(1,length(a)),covm,N){
	snp_row<-1:length(beta_a[,1])
	
	i=1
	for (i in 1:length(var_y)){
		beta_a[,i]=beta_a[,i]/sqrt(var_y[i])
		se[,i]=se[,i]/sqrt(var_y[i])
	}
	
	covm=(1/(sqrt(diag(covm))%o%sqrt(diag(covm))))*covm
	
	N_min_index=sapply(snp_row, function(x) which.min(N[x,]))
	se_N_min<-sapply(snp_row, function (x) se[x,N_min_index[x]])
	N_min<-sapply(snp_row, function(x) N[x,N_min_index[x]])
	beta_N_min<-sapply(snp_row, function(x) beta_a[x,N_min_index[x]])
	
	#var_y_N_min<-sapply(snp_row, function(x) var_y[N_min_index[x]])
	
	var_y_a=sum(covm*(a%o%a))
	b=(beta_a%*%a)
	
	var_y_N_min_to_var_g_N_min_ratio=se_N_min^2+(beta_N_min^2)/N_min
	
	#se2=var_y_N_min_to_var_g_N_min_ratio*(var_y_a/var_y_N_min)
	se2=var_y_N_min_to_var_g_N_min_ratio*(var_y_a/1)
	
	se2=se2-b^2/N_min
	
	se=sqrt(se2)
	
	b=b/sqrt(var_y_a)
	se=se/sqrt(var_y_a)
	
	out=cbind(b=b,se=se)
	colnames(out)=c("b","se")
	out=as.data.frame(out)
	return(out)
}


sh_gwas=GWAS_linear_combination(a=as.numeric(aa[2,]),beta_a=betas,se=se,covm=as.matrix(covm),N=sample_size)

sh_gwas=mutate(x,Z=b/se,p=pchisq(Z^2,1,low=F))
sh_gwas=mutate(x,SNP=gwas[[1]]$rs_id)
sh_gwas=mutate(x,A1=gwas[[1]]$ea,A2=gwas[[1]]$ra,N=336107,chr=gwas[[1]]$chr,pos=gwas[[1]]$bp,
			 eaf=gwas[[1]]$eaf)

head(sh_gwas,n=2)

data.table::fwrite(
	x, 
	row.names=F,
	file = '/../data/anthropometry_results/four_traits/linear_combination/SH_GWAS.txt',
	sep = '\t')
  }
  
