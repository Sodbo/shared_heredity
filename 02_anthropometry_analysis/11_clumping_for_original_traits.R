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
	lt <- function_for_shlop_28_12_2017(sst_sm,p_value="p",pos="bp",snp="rs_id", delta=5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
		lt=cbind(lt,trait)
		out=rbind(out,lt)
	}
}
	#Clumping for SH
	
	sst_file <-'../../data/01_anthropometry_results/linear_combination/SH_GWAS.txt'
	sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
	sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
	trait <- 'SH'
	lt <- function_for_shlop_28_12_2017(sst_sm,p_value="p",pos="pos",snp="SNP", delta=5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
		#correct names to merge tables
		colnames(lt)=c('beta','se','z','p','rs_id','ea','ra','n','chr','bp','eaf')
		lt=cbind(lt,trait)
		out<-out[c(colnames(lt))]
		out=rbind(out,lt)
		}
	
dim(out)

bt <- function_for_shlop_28_12_2017(out,trait="trait",p_value="p",pos="bp",snp="rs_id",delta=5e5,chr="chr")
bt <- bt[order(bt$chr,bt$bp),]
colnames(bt)
str(bt)
bt$trait <- as.character(bt$trait)
str(bt)
write.csv(bt, paste0(path, result_file_name))

