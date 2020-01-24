# This script fuse tables obtained by grep of clumped SNP for Shared heredity gwas.
path='../data/anthropometry_results/four_traits/GWAS/scaled_filtered/'
trait4049<-read.csv(paste0(path,'ID_4049_clump.txt'),header=F, stringsAsFactors=F)
trait4050<-read.csv(paste0(path,'ID_4050_clump.txt'),header=F, stringsAsFactors=F)
trait4058<-read.csv(paste0(path,'ID_4058_clump.txt'),header=F, stringsAsFactors=F)
trait4179<-read.csv(paste0(path,'ID_4179_clump.txt'),header=F, stringsAsFactors=F)
sh<-read.delim(paste0(path,'SH_clump.txt'),header=F, stringsAsFactors=F)
clump<-read.csv(paste0(path,'Clumping_for_all_orignal_traits_and_SH.csv'))

#remove empty rs_id
trait4049<-trait4049[-122,]


clump$p_4049<-trait4049[,12]
clump$p_4050<-trait4050[,12]
clump$p_4058<-trait4058[,12]
clump$p_4179<-trait4179[,12]
clump$p_SH<-sh[,4]
clump$se_4049<-trait4049[,11]
clump$se_4050<-trait4050[,11]
clump$se_4058<-trait4058[,11]
clump$se_4179<-trait4179[,11]
clump$se_SH<-sh[,2]
clump$beta_4049<-trait4049[,10]
clump$beta_4050<-trait4050[,10]
clump$beta_4058<-trait4058[,10]
clump$beta_4179<-trait4179[,10]
clump$beta_SH<-sh[,1]

write.csv(clump,paste0(path,'summary_for_SH_clumping.csv'))
