# Aim of this script is to do clumping for original traits and SH and generate a final table

library(data.table)

setwd('/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases/')
source("../00_core_functions/clumping.R", chdir = F)

path <- "../../../data/03_neurodegenerative_diseases/several_traits/"
input_file_name<-c('../../../data/03_neurodegenerative_diseases/BIP/03_gc_corrected/bip_gc_corrected.csv',
		   '../../../data/03_neurodegenerative_diseases/MDD/03_gc_corrected/mdd_gc_corrected.csv',
		   '../../../data/03_neurodegenerative_diseases/SCZ/03_gc_corrected/scz_gc_corrected.csv',
		   '/../../../data/03_neurodegenerative_diseases/several_traits/linear_combination/03_gc_corrected/sh_gc_corrected.csv')
result_file_name<-'Clumping_for_all_orignal_traits_and_SH.txt'
out <- NULL
thr <- 5e-8


for (input in input_file_name) {
	sst_file <- paste0(path,input)
	sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
	sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
	trait<-strsplit(input, '03_gc_corrected/')[[1]][2]
	trait <- strsplit(trait, '_')[[1]][1]
	lt <- function_for_shlop_28_12_2017(sst_sm,p_value="p_gc",pos="bp",snp="rs_id", delta=5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
		lt=cbind(lt,trait)
		out=rbind(out,lt)
	}
}

dim(out)

bt <- function_for_shlop_28_12_2017(out,trait="trait",p_value="p_gc",pos="bp",snp="rs_id",delta=5e5,chr="chr")
bt <- bt[order(bt$chr,bt$bp),]
colnames(bt)
str(bt)
bt$trait <- as.character(bt$trait)
str(bt)
write.csv(bt, paste0(path, result_file_name))

