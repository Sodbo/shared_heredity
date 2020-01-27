library(data.table)
library(dplyr)

path_to_result_directory<-'../../data/anthropometry_results/four_traits/GWAS/scaled_filtered/'
gwas_files<-list.files(path_to_result_directory, full.names=T, pattern="ID_\\d+.csv")
gwas<-lapply(gwas_files, fread)

aa<-read.table('../../data/anthropometry_results/four_traits/alphas.txt', row.names=1)
covm<-read.table('../../data/anthropometry_results/four_traits/pheno_corr_matrix.txt', row.names=1,  check.names=F)

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
	
	out=cbind(b=b,se=se,N=N_min)
	colnames(out)=c("b","se","N")
	out=as.data.frame(out)
	return(out)
}

cov_gi_alpha=function(a,i,covm=covm){
  cov_gi_sum_giai=sum(covm[i,]*a)
  var_gi=covm[i,i]
  var_sum_giai=sum(covm*(a%o%a)) #gen(a)
  
  cor_g_a=cov_gi_sum_giai
  return(cor_g_a)
} 

alphas=as.numeric(aa[2,])
b1=cov_gi_alpha(alphas,1,as.matrix(covm))

tr1_sh_gwas=GWAS_linear_combination(a=c(1,0,0,0)-alphas*b1,beta_a=betas,se=se,covm=as.matrix(covm),N=sample_size)

tr1_sh_gwas=mutate(tr1_sh_gwas,Z=b/se,p=pchisq(Z^2,1,low=F))
tr1_sh_gwas=mutate(tr1_sh_gwas,SNP=gwas_reordered[[1]]$rs_id)
tr1_sh_gwas=mutate(tr1_sh_gwas,A1=gwas_reordered[[1]]$ea,A2=gwas_reordered[[1]]$ra,chr=gwas_reordered[[1]]$chr,pos=gwas_reordered[[1]]$bp,
			 eaf=gwas_reordered[[1]]$eaf)

head(tr1_sh_gwas,n=2)
data.table::fwrite(
	tr1_sh_gwas, 
	row.names=F,
	file = '../../data/anthropometry_results/four_traits/linear_combination/Tr4049-SH_GWAS.txt',
	sep = '\t')
	
b2=cov_gi_alpha(alphas,2,as.matrix(covm))

tr2_sh_gwas=GWAS_linear_combination(a=c(0,1,0,0)-alphas*b2,beta_a=betas,se=se,covm=as.matrix(covm),N=sample_size)

tr2_sh_gwas=mutate(tr2_sh_gwas,Z=b/se,p=pchisq(Z^2,1,low=F))
tr2_sh_gwas=mutate(tr2_sh_gwas,SNP=gwas_reordered[[1]]$rs_id)
tr2_sh_gwas=mutate(tr2_sh_gwas,A1=gwas_reordered[[1]]$ea,A2=gwas_reordered[[1]]$ra,chr=gwas_reordered[[1]]$chr,pos=gwas_reordered[[1]]$bp,
			 eaf=gwas_reordered[[1]]$eaf)

head(tr2_sh_gwas,n=2)
data.table::fwrite(
	tr2_sh_gwas, 
	row.names=F,
	file = '../../data/anthropometry_results/four_traits/linear_combination/Tr4050-SH_GWAS.txt',
	sep = '\t')

b3=cov_gi_alpha(alphas,3,as.matrix(covm))

tr3_sh_gwas=GWAS_linear_combination(a=c(0,0,1,0)-alphas*b3,beta_a=betas,se=se,covm=as.matrix(covm),N=sample_size)

tr3_sh_gwas=mutate(tr3_sh_gwas,Z=b/se,p=pchisq(Z^2,1,low=F))
tr3_sh_gwas=mutate(tr3_sh_gwas,SNP=gwas_reordered[[1]]$rs_id)
tr3_sh_gwas=mutate(tr3_sh_gwas,A1=gwas_reordered[[1]]$ea,A2=gwas_reordered[[1]]$ra,chr=gwas_reordered[[1]]$chr,pos=gwas_reordered[[1]]$bp,
			 eaf=gwas_reordered[[1]]$eaf)

head(tr3_sh_gwas,n=2)
data.table::fwrite(
	tr3_sh_gwas, 
	row.names=F,
	file = '../../data/anthropometry_results/four_traits/linear_combination/Tr4058-SH_GWAS.txt',
	sep = '\t')

b4=cov_gi_alpha(alphas,4,as.matrix(covm))

tr4_sh_gwas=GWAS_linear_combination(a=c(0,0,0,1)-alphas*b4,beta_a=betas,se=se,covm=as.matrix(covm),N=sample_size)

tr4_sh_gwas=mutate(tr4_sh_gwas,Z=b/se,p=pchisq(Z^2,1,low=F))
tr4_sh_gwas=mutate(tr4_sh_gwas,SNP=gwas_reordered[[1]]$rs_id)
tr4_sh_gwas=mutate(tr4_sh_gwas,A1=gwas_reordered[[1]]$ea,A2=gwas_reordered[[1]]$ra,chr=gwas_reordered[[1]]$chr,pos=gwas_reordered[[1]]$bp,
			 eaf=gwas_reordered[[1]]$eaf)

head(tr4_sh_gwas,n=2)
data.table::fwrite(
	tr4_sh_gwas, 
	row.names=F,
	file = '../../data/anthropometry_results/four_traits/linear_combination/Tr4179-SH_GWAS.txt',
	sep = '\t')
