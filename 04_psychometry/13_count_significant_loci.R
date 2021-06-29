# This script parses joint clumping results to extract number of significant loci for each trait
trait_designation<-c('bip','mdd','scz','happiness','SH')
library(data.table)
clump<-fread('../../data/03_neurodegenerative_diseases/several_traits/four_traits/clumping_of_SH_and_orig_traits_partIV_5e-8.txt')
significant<-lapply(trait_designation, function (x) grepl(paste0('\\b',x,'\\b'),clump$traits, perl=T))
count<-sapply(significant, table)
colnames(count)<-trait_designation
count[2,] 
 
