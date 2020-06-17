# This script parses joint clumping results to extract number of significant loci for each trait
trait_designation<-c(1287001,1287003,1287004,'SH')
library(data.table)
clump<-fread('../../data/02_Lipids/three_traits/clumping_of_SH_and_orig_traits_partIV_5e-8.txt')
significant<-lapply(trait_designation, function (x) grepl(x,clump$traits))
count<-sapply(significant, table)
colnames(count)<-trait_designation
count[2,] 
