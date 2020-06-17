# This script parses joint clumping results to extract number of significant loci for each trait
trait_designation<-c(191,192,193,194,199,'SH')
library(data.table)
clump<-fread('../../data/01_anthropometry_results/five_traits/joint_clumping_and_enrichment_new/clumping_of_SH_and_orig_traits_partIV_5e-8.txt')
significant<-lapply(trait_designation, function (x) grepl(x,clump$traits))
count<-sapply(significant, table)
colnames(count)<-trait_designation
count[2,]

