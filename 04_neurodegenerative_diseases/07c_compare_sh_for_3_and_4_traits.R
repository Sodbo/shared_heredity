# Aim of this script is to compare SH obtained for 3 ttraits with 
# the new one for four traits

library(data.table)
library(dplyr)

path <- '../../../data/03_neurodegenerative_diseases/several_traits/four_traits/'

sh <- fread(paste0(path,'linear_combination/00_raw_data/SH_GWAS.txt'))
sh_old <- fread('../../../data/03_neurodegenerative_diseases/several_traits/linear_combination/00_raw_data/SH_GWAS.txt')

table(sh$SNP %in% sh_old$SNP)
ind <- match(sh$SNP, sh_old$SNP)
table(sh$SNP == sh_old$SNP[ind])
sh_old <- sh_old[ind, ]

png(paste0(path, 'beta.png'))
plot(sh$b[1:10000], sh_old$b[1:10000], xlab='new', ylab='old', main='beta')
dev.off()
png(paste0(path,'se.png'))
plot(sh$se[1:10000], sh_old$se[1:10000], xlab='new', ylab='old', main='se')
dev.off()
png(paste0(path,'Z.png'))
plot(sh$Z[1:10000], sh_old$Z[1:10000], xlab='new', ylab='old', main='Z-score')
dev.off()

