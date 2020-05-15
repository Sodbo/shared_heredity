# Aim of this script is to obtain GWAS summary statistics for shared heredity

library('data.table')
library(dplyr)
source("../00_core_functions/linear_combination_v3.R")

path_to_gwas_directory<-'../../data/01_anthropometry_results/five_traits/GWAS/'
gwas_files<-list.files(path_to_gwas_directory, full.names=T, pattern="ID_\\d+.csv")
gwas<-lapply(gwas_files, fread)

gip1<-read.table('../../data/01_anthropometry_results/five_traits/GIP1.txt', row.names=1)
covm<-read.table('../../data/01_anthropometry_results/five_traits/pheno_corr_matrix.txt', row.names=1,  check.names=F)

rs_id<-lapply(gwas, function(x) x$rs_id)
snps<-Reduce(intersect,rs_id)
ind<-lapply(rs_id, function(x) match(snps,x))

#Is it necessary to remove NA for the case of not whole intersection?
#ind<-lapply(ind,function(x) )

gwas_reordered<-lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]],])
z <- sapply(gwas_reordered, function(x) x$z)
sample_size <- sapply(gwas_reordered, function(x) x$n)

#eaf is better to use for gwas with max sample size, in this case sample sizes 

gip1_gwas<-GWAS_linear_combination_Z_based(a=as.numeric(gip1[,1]),Z = z,covm=as.matrix(covm),N=sample_size, eaf=gwas_reordered[[1]]$eaf)

gip1_gwas=mutate(gip1_gwas,Z=b/se,p=pchisq(Z^2,1,low=F))
gip1_gwas=mutate(gip1_gwas,SNP=gwas_reordered[[1]]$rs_id)
gip1_gwas=mutate(gip1_gwas,A1=gwas_reordered[[1]]$ea,A2=gwas_reordered[[1]]$ra,chr=gwas_reordered[[1]]$chr,pos=gwas_reordered[[1]]$bp,
			 eaf=gwas_reordered[[1]]$eaf)

head(gip1_gwas,n=2)
dir.create('../../data/01_anthropometry_results/five_traits/GIP1/')
data.table::fwrite(
	gip1_gwas, 
	row.names=F,
	file = '../../data/01_anthropometry_results/five_traits/GIP1/GIP1_GWAS.txt',
	sep = '\t')
