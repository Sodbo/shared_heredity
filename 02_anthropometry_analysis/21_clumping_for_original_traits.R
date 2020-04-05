# This script is for joint clumping for the original traits without
# shared heredity, to find shared hits - SNPs showing significance
# at low threshold,
thr <- 1e-5

library(data.table) 
source('../00_core_functions/clumping.R')
path <- "../../data/01_anthropometry_results/GWAS/scaled_filtered/"
input_file_name<-c('ID_4049_gc_corrected.csv','ID_4050_gc_corrected.csv','ID_4058_gc_corrected.csv','ID_4179_gc_corrected.csv')
result_file_name<-'shared_hits/Clumping_for_all_orignal_traits.txt'
out <- NULL


for (input in input_file_name) {
	sst_file <- paste0(path,input)
	sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
	sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
	trait<-input
	lt <- function_for_shlop_29_03_2020(sst_sm,p_value="p_gc",pos="bp",snp="rs_id", delta=5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
		lt=cbind(lt,trait)
		out=rbind(out,lt)
	}
}

bt <- function_for_shlop_29_03_2020(out,trait="trait",p_value="p_gc",pos="bp",snp="rs_id",delta=5e5,chr="chr")
bt <- bt[order(bt$chr,bt$bp),]
colnames(bt)
str(bt)
bt$trait <- as.character(bt$trait)
str(bt)
write.csv(bt, paste0(path, result_file_name))
fwrite(as.list(bt$rs_id), paste0(path, 'shared_hits/list_of_clumped_rs_id.txt'),sep='\n')
