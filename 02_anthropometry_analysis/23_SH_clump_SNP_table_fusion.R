# This script add column p_val_SH containing shared heredity p-value,
# draw boxplot for SH p-value against shared hit level
thr=5e-8
path<-'../../data/01_anthropometry_results/GWAS/scaled_filtered/shared_hits/'
sh<-read.csv(paste0(path,'SH_clumped_snps_v2.txt'),header=F, stringsAsFactors=F)
clump<-read.csv(paste0(path,'Clumping_for_all_orignal_traits.txt'),stringsAsFactors=F)
# read file to get column_names for sum stats
sumstat_header_sh<-read.csv('../../data/01_anthropometry_results/linear_combination/02_unification_out/SH/SH_done.csv', nrows=1)

colnames(sh)<-colnames(sumstat_header_sh)

#grepping of shared hits is not in right order, so it necessary to order it.
ind <- match(clump$rs_id, sh$rs_id)
sh<-sh[ind,]

#Add coulumns with SH p-value and with 
clump$p_val_SH<-sh$p
clump$SH_significan<-sh$p<thr
write.csv(clump,paste0(path,'summary_for_shared_hits_and_SH.csv'))

png(paste0(path,'Significance_of_SH_for_shared_hits.png'))
boxplot(-log10(p_val_SH) ~ Ntraits, data = clump, ylim=c(0,40) , xlab = 'N of traits', ylab = '-log_10(p-value) for SH')
dev.off()

reg<-lm(-log(p_val_SH)~Ntraits, data=clump)
corr<-cor(clump$Ntraits, -log10(clump$p_val_SH))
det=corr^2
cortest<-cor.test(clump$Ntraits, -log10(clump$p_val_SH), method='kendal' )

sink(paste0(path,"regression.txt"))
print(summary(reg))
print(paste0('r=',corr))
print(paste0('r^2=',det))
print(cortest)
sink()
