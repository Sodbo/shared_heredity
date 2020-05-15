#checking if old results is differnt from new ones.
#
path<-'../../data/01_anthropometry_results/five_traits/'
library(data.table)

sh<-fread(paste0(path,'linear_combination/SH_GWAS.txt'))
sh_old<-fread('../../data/01_anthropometry_results/linear_combination/SH_GWAS.txt')
png(paste0(path,'linear_combination/validity_check/beta.png'))
plot(sh$b, sh_old$b, xlab='new', ylab='old', main='beta')
dev.off()
png(paste0(path,'linear_combination/validity_check/se.png'))
plot(sh$se, sh_old$se, xlab='new', ylab='old', main='se')
dev.off()
png(paste0(path,'linear_combination/validity_check/Z.png'))
plot(sh$Z,sh_old$Z)
dev.off()

source('../00_core_functions/heritability_of_linear_combination.R')
covm<-read.delim(paste0(path,'gene_cov_matrix.txt'),sep=' ', row.names=1)
phem<-read.delim(paste0(path,'pheno_corr_matrix.txt'),sep=' ', row.names=1)
a<-read.delim(paste0(path,'alphas.txt'),sep='', row.names=1)
a<-as.numeric(a[2,])

H2(a,covm=covm,phem=phem)
