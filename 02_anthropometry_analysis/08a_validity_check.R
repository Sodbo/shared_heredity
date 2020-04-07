#checking if old results is differnt from new ones.
#

library(data.table)

sh<-fread('../../data/01_anthropometry_results/linear_combination/SH_GWAS.txt')
sh_old<-fread('../../data/01_anthropometry_results/linear_combination/old_results/SH_GWAS.txt')
png('../../data/01_anthropometry_results/linear_combination/old_results/new_results_validity/beta.png')
plot(sh$b, sh_old$b, xlab='new', ylab='old', main='beta')
dev.off()
png('../../data/01_anthropometry_results/linear_combination/old_results/new_results_validity/se.png')
plot(sh$se, sh_old$se, xlab='new', ylab='old', main='se')
dev.off()
png('../../data/01_anthropometry_results/linear_combination/old_results/new_results_validity/Z.png')
plot(sh$Z,sh_old$Z)
dev.off()

source('../00_core_functions/heritability_of_linear_combination.R')
covm<-read.delim('../../data/01_anthropometry_results/gene_cov_matrix.txt',sep=' ', row.names=1)
phem<-read.delim('../../data/01_anthropometry_results/pheno_corr_matrix.txt',sep=' ', row.names=1)
a<-read.delim('../../data/01_anthropometry_results/alphas.txt',sep=' ', row.names=1)
a<-as.numeric(a[2,1:4])

H2(a,covm=covm,phem=phem)
