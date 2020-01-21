# This script fuse tables obtained by grep of clumped SNP for Shared heredity gwas.
path='../data/anthropometry_results/four_traits/GWAS/scaled_filtered/'
g4049<-read.csv(paste0(path,'ID_4049_SH_clump.txt'),header=F, stringsAsFactors=F)
g4050<-read.csv(paste0(path,'ID_4050_SH_clump.txt'),header=F, stringsAsFactors=F)
g4058<-read.csv(paste0(path,'ID_4058_SH_clump.txt'),header=F, stringsAsFactors=F)
g4179<-read.csv(paste0(path,'ID_4179_SH_clump.txt'),header=F, stringsAsFactors=F)
SH<-read.csv(paste0(path,'../../linear_combination/clumping_results_1000kb.txt'))



SH$p_4049<-g4049[,12]
SH$p_4050<-g4050[,12]
SH$p_4058<-g4058[,12]
SH$p_4179<-g4179[,12]

write.csv(SH,paste0(path,'summary_for_SH_clumping.csv'))
