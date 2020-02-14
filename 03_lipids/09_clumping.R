#Clumping for the SH GWAS

library(data.table) 
source("../00_core_functions/clumping.R")

path <- "/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/linear_combination/"
out <- NULL
thr <- 5e-8


sst_file <- paste0(path,'SH_Lipids_GWAS.txt')
sst <- fread(sst_file, header=T, stringsAsFactors = F, data.table=F)
sst_sm <- sst[(pmin(sst$eaf,1-sst$eaf)>=0.01),]
trait <- 'SH for anthropometric traits'
lt <- function_for_shlop_28_12_2017(sst_sm,p_value="p",pos="pos",snp="SNP", delta=5e5,chr="chr", thr=thr)
if (nrow(lt)>0 ) {
		lt=cbind(lt,trait)
		out=rbind(out,lt)
	}

dim(out)

write.csv(out, paste0(path,'clumping_results_sh_lipids.csv'))
