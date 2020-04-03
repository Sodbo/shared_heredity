# This script fuse tables obtained by grep of clumped SNP for Shared heredity gwas.
path='../../data/01_anthropometry_results/GWAS/scaled_filtered/'
trait4049<-read.delim(paste0(path,'ID_4049_clumped_with_SH.txt'),header=F, stringsAsFactors=F)
trait4050<-read.delim(paste0(path,'ID_4050_clumped_with_SH.txt'),header=F, stringsAsFactors=F)
trait4058<-read.delim(paste0(path,'ID_4058_clumped_with_SH.txt'),header=F, stringsAsFactors=F)
trait4179<-read.delim(paste0(path,'ID_4179_clumped_with_SH.txt'),header=F, stringsAsFactors=F)
sh<-read.delim(paste0(path,'SH_clumped_with_traits.txt'),header=F, stringsAsFactors=F)
clump<-read.csv(paste0(path,'Clumping_for_all_orignal_traits_and_SH.txt'),stringsAsFactors=F)
# read file to get column_names for sum stats
sumstat_header<-read.delim(paste0(path,'ID_4049_gc_corrected.csv'), nrows=1)
sumstat_header_sh<-read.delim('../../data/01_anthropometry_results/linear_combination/SH_GWAS.txt', nrows=1)

colnames(trait4049)<-colnames(sumstat_header)
colnames(trait4050)<-colnames(sumstat_header)
colnames(trait4058)<-colnames(sumstat_header)
colnames(trait4179)<-colnames(sumstat_header)
colnames(sh)<-colnames(sumstat_header_sh)

#After grepping with previous script (14) the extra line repeating the first line appeared, the reason is not know yet (probably due to >> operator). And that's the following lines of code is for deleting last line in all the clumped files:
trait4049<-trait4049[1:684,]
trait4050<-trait4050[1:684,]
trait4058<-trait4058[1:684,]
trait4179<-trait4179[1:684,]
sh<-sh[39:(39+683),]

clump$p_4049<-trait4049$p_gc
clump$p_4050<-trait4050$p_gc
clump$p_4058<-trait4058$p_gc
clump$p_4179<-trait4179$p_gc
if(is.null(sh$p_gc)){
	sh$p_gc<-sh$p
}
clump$p_SH<-sh$p_gc

clump$se_4049<-trait4049$se
clump$se_4050<-trait4050$se
clump$se_4058<-trait4058$se
clump$se_4179<-trait4179$se
clump$se_SH<-sh$se

clump$beta_4049<-trait4049$beta
clump$beta_4050<-trait4050$beta
clump$beta_4058<-trait4058$beta
clump$beta_4179<-trait4179$beta
clump$beta_SH<-sh$b

write.csv(clump,paste0(path,'summary_for_SH_clumping.csv'))
