#checking if old results is differnt from new ones.
#
path<-'../../data/01_anthropometry_results/five_traits/'
library(data.table)

sh<-fread(paste0(path,'linear_combination/SH_GWAS.txt'))
maxH<-fread(paste0(path,'linear_combination/maxH_GWAS.txt'))
png(paste0(path,'linear_combination/validity_check_maxH/beta.png'))
plot(sh$b, maxH$b, xlab='sh', ylab='maxH', main='beta')
abline(0,1)
dev.off()
png(paste0(path,'linear_combination/validity_check_maxH/se.png'))
plot(sh$se, maxH$se, xlab='sh', ylab='maxH', main='se')
abline(0,1)
dev.off()
png(paste0(path,'linear_combination/validity_check_maxH/Z.png'))
plot(sh$Z,maxH$Z, xlab='sh', ylab='maxH', main='Z')
abline(0,1)
dev.off()

source('../00_core_functions/heritability_of_linear_combination.R')
covm<-read.delim(paste0(path,'gene_cov_matrix.txt'),sep=' ', row.names=1)
phem<-read.delim(paste0(path,'pheno_corr_matrix.txt'),sep=' ', row.names=1)
a<-read.delim(paste0(path,'alphas.txt'),sep='', row.names=1)
a<-as.numeric(a[,1])

H2(a,covm=covm,phem=phem)
