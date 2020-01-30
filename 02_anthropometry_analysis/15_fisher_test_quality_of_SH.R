clump<-read.csv('../../data/anthropometry_results/four_traits/GWAS/scaled_filtered/summary_for_SH_clumping.csv', stringsAsFactors=F) 
#Add new zero columns
clump$all_significant_and_SH=0
clump$not_significant_and_SH=0 
clump$all_significant_not_SH=0
clump$not_significant_not_SH=0 


clump$is_SH<-grepl('\\bSH\\b',clump$traits,perl=T) #SH is significant

clump$all_significant_and_SH[clump$Ntraits==5]<-1 #all traits are significant including SH
clump$not_significant_and_SH[clump$Ntraits!=5 & clump$is_SH==T]<-1 #not all traits are significant but SH is significant
clump$all_significant_not_SH[clump$Ntraits==4 & clump$is_SH==F]<-1 #all traits are significant excluding SH
clump$not_significant_not_SH[clump$Ntraits!=5 & clump$is_SH==F]<-1 #not all significant and SH is not significant
n11=sum(clump$all_significant_and_SH)
n12=sum(clump$not_significant_and_SH)
n21=sum(clump$all_significant_not_SH)
n22=sum(clump$not_significant_not_SH)
contigency_table=matrix(c(n11,n12,n21,n22),nrow=2,ncol=2,byrow=T)
colnames(contigency_table)<-c('all_significant','not all significant')
rownames(contigency_table)<-c('SH is significant','SH is not significant')
write.table(contigency_table, '../../data/anthropometry_results/four_traits/GWAS/scaled_filtered/contigency_table.csv')

result<-fisher.test(contigency_table)

write.table(result$p.value,'../../data/anthropometry_results/four_traits/GWAS/scaled_filtered/SH_enrichment_analysis.txt') 
')
