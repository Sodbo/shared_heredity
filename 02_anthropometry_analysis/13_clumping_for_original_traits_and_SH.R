#Clumping for the original traits

library(data.table) 
source('../00_core_functions/clumping.R')
path <- "../../data/01_anthropometry_results/GWAS/scaled_filtered/"
input_file_name<-c('ID_4049_gc_corrected.csv','ID_4050_gc_corrected.csv','ID_4058_gc_corrected.csv','ID_4179_gc_corrected.csv')
result_file_name<-'Clumping_for_all_orignal_traits_and_SH.txt'
out <- NULL
thr <- 5e-8

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
	#Clumping for SH
	
	sst_file <-'../../data/01_anthropometry_results/linear_combination/02_unification_out/SH/SH_done.csv'
	sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
	sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
	if(is.null(sst_sm$p_gc)){ # if genomic control is not done
		sst_sm$p_gc<-sst_sm$p # add column to 'p_gc' to correctly declare p_value column in the followng joint clumping process of the SH and original traits
	}
	lt <- function_for_shlop_29_03_2020(sst_sm,p_value="p_gc",pos="bp",snp="rs_id", delta=5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
		lt$trait <- 'SH'
		out<-out[c(colnames(lt))]
		out=rbind(out,lt)
	}

dim(out)

bt <- function_for_shlop_29_03_2020(out,trait="trait",p_value="p_gc",pos="bp",snp="rs_id",delta=5e5,chr="chr")
bt <- bt[order(bt$chr,bt$bp),]
colnames(bt)
str(bt)
bt$trait <- as.character(bt$trait)
str(bt)
write.csv(bt, paste0(path, result_file_name))
fwrite(as.list(bt$rs_id), paste0(path, 'list_of_clumped_rs_id.txt'),sep='\n')
