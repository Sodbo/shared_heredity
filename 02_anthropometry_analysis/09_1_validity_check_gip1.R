#checking if old results is differnt from new ones.
#
path<-'../../data/01_anthropometry_results/five_traits/'
library(data.table)

sh<-fread(paste0(path,'linear_combination/SH_GWAS.txt'))
gip1<-fread(paste0(path,'GIP1/GIP1_GWAS.txt'))
png(paste0(path,'GIP1/validity_check/beta.png'))
plot(sh$b, gip1$b, xlab='sh', ylab='gip1', main='beta')
abline(0,1)
dev.off()
png(paste0(path,'GIP1/validity_check/se.png'))
plot(sh$se, gip1$se, xlab='sh', ylab='gip1', main='se')
abline(0,1)
dev.off()
png(paste0(path,'GIP1/validity_check/Z.png'))
plot(sh$Z,gip1$Z, xlab='sh', ylab='gip1', main='Z')
abline(0,1)
dev.off()

source('../00_core_functions/heritability_of_linear_combination.R')
covm<-read.delim(paste0(path,'gene_cov_matrix.txt'),sep=' ', row.names=1)
phem<-read.delim(paste0(path,'pheno_corr_matrix.txt'),sep=' ', row.names=1)
a<-read.delim(paste0(path,'GIP1.txt'),sep='', row.names=1)
a<-as.numeric(a[,1])

H2(a,covm=covm,phem=phem)
