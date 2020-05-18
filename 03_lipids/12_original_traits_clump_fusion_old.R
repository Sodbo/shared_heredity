#Fusion of clumping results 

library(data.table)

path='../data/anthropometry_results/four_traits/GWAS/scaled_filtered/' 

clump_files<-list.files(path, full.names=T, pattern="Clumping_for_\\d+.txt")
clump<-lapply(clump_files, fread)

rs_id<-sapply(clump, function(x) x[,1])
rs<-unique(unlist(rs_id))
#tr4049<-read.csv(paste0(path,'Clumping_for_4049.txt'), stringsAsFactors=F)
#tr4050<-read.csv(paste0(path,'Clumping_for_4050.txt'), stringsAsFactors=F)
#tr4058<-read.csv(paste0(path,'Clumping_for_4058.txt'), stringsAsFactors=F)
#tr4179<-read.csv(paste0(path,'Clumping_for_4179.txt'), stringsAsFactors=F)
SH<-read.csv(paste0(path,'../../linear_combination/clumping_results_1000kb.txt'))

union(tr4049[,1],)
