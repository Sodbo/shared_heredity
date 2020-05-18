# This script fuse tables obtained by grep of clumped SNP for Shared heredity gwas.
path='/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/GWAS/scaled_filtered/'
g1287001<-read.csv(paste0(path,'ID_1287001_SH_clump.txt'),header=F, stringsAsFactors=F)
g1287003<-read.csv(paste0(path,'ID_1287003_SH_clump.txt'),header=F, stringsAsFactors=F)
g1287004<-read.csv(paste0(path,'ID_1287004_SH_clump.txt'),header=F, stringsAsFactors=F)
sh<-read.delim(paste0(path,'SH_clump.txt'),header=F, stringsAsFactors=F)
clump <-read.csv(paste0(path,'Clumping_for_3_orignal_traits_and_SH.csv'))

clump$p_1287001<-g1287001[,12]
clump$p_1287003<-g1287003[,12]
clump$p_1287004<-g1287004[,12]
clump$p_SH<-sh[,5]

clump$se_1287001<-g1287001[,11]
clump$se_1287003<-g1287003[,11]
clump$se_1287004<-g1287004[,11]
clump$se_SH<-sh[,2]

clump$beta_1287001<-g1287001[,10]
clump$beta_1287003<-g1287003[,10]
clump$beta_1287004<-g1287004[,10]
clump$beta_SH<-sh[,1]


write.csv(clump,paste0(path,'summary_for_SH_clumping.csv'))
