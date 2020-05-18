#Clumping for the original traits

library(data.table) 

path <- "/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/GWAS/scaled_filtered/"
input_file_name<-c('ID_1287001.csv','ID_1287003.csv','ID_1287004.csv')
result_file_name<-'Clumping_for_3_orignal_traits_and_SH.csv'
out <- NULL
thr <- 5e-8
source("../00_core_functions/clumping.R")


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
	
	sst_file <-'/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/linear_combination/SH_Lipids_GWAS_2.txt'
	sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
	sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
	trait <- 'SH'
	lt <- function_for_shlop_28_12_2017(sst_sm,p_value="p",pos="pos",snp="SNP", delta=5e5,chr="chr", thr=thr)
	if (nrow(lt)>0 ) {
		#correct names to merge tables
		colnames(lt)=c('beta','se','n','z','p','rs_id','ea','ra','chr','bp','eaf')
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

write.csv(
	bt, 
	file = paste0(path, result_file_name))


