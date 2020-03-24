library(data.table)
library(dplyr)
source("../00_core_functions/linear_combination.r")

gwas_files<-c('../../../data/neurodegenerative_diseases/BIP/02_unification_results/bip_output_done.csv',
	      '../../../data/neurodegenerative_diseases/MDD/02_unification_results/mdd_output_done.csv',
	      '../../../data/neurodegenerative_diseases/SCZ/02_unification_results/scz_output_done.csv')
gwas<-lapply(gwas_files, fread)

aa<-read.table('../../../data/neurodegenerative_diseases/several_traits/alphas.txt', row.names=1)
covm<-read.table('../../../data/neurodegenerative_diseases/several_traits/pheno_corr_matrix.txt', row.names=1,  check.names=F)

rs_id<-lapply(gwas, function(x) x$rs_id)
snps<-Reduce(intersect,rs_id)
ind<-lapply(rs_id, function(x) match(snps,x))

#Is it necessary to remove NA for the case of not whole intersection?
#ind<-lapply(ind,function(x) )

gwas_reordered<-lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]],])
betas<-sapply(gwas_reordered, function(x) x$beta)
se <- sapply(gwas_reordered, function(x) x$se)
sample_size <- sapply(gwas_reordered, function(x) x$n)

sh_gwas=GWAS_linear_combination(a=as.numeric(aa[2,]),beta_a=betas,se=se,covm=as.matrix(covm),N=sample_size)

sh_gwas=mutate(sh_gwas,Z=b/se,p=pchisq(Z^2,1,low=F))
sh_gwas=mutate(sh_gwas,SNP=gwas_reordered[[1]]$rs_id)
sh_gwas=mutate(sh_gwas,A1=gwas_reordered[[1]]$ea,A2=gwas_reordered[[1]]$ra,chr=gwas_reordered[[1]]$chr,pos=gwas_reordered[[1]]$bp,
			 eaf=gwas_reordered[[1]]$eaf)

head(sh_gwas,n=2)
dir.create('../../../data/neurodegenerative_diseases/several_traits/linear_combination/')
data.table::fwrite(
	sh_gwas, 
	row.names=F,
	file = '../../../data/neurodegenerative_diseases/several_traits/linear_combination/SH_GWAS.txt',
	sep = '\t')
