#This script makes new summary statistics file with column p_gc containing p-values corrected for genomic control.
# To run the script be sure the order of files in directory correspondes the order of gwas_ids in .csv files

sh_gwas_id=201

library(data.table)
path_to_result_directory='../../data/01_anthropometry_results/five_traits/GWAS/'
gwas_files<-list.files(path_to_result_directory, full.names=T, pattern="ID_\\d+.csv")
intercept_data<-read.csv(paste0(path_to_result_directory,'intrecept_data.csv')) 
gwas<-lapply(gwas_files, fread)

sh<-fread('../../data/01_anthropometry_results/five_traits/linear_combination/SH_GWAS.txt')
gwas<-c(gwas,list(sh))

names(gwas)<-regmatches(gwas_files, regexpr('(?<=ID_)(\\d+)(?=\\.csv)', gwas_files, perl=T))
names(gwas)[length(gwas)]=sh_gwas_id

chi=lapply(gwas, function(x) qchisq(x$p, df=1, lower.tail=F))
chi_gc=lapply(1:length(chi), function(x) chi[[x]]/intercept_data$intercept[intercept_data$gwas_id==names(gwas)[x]])
gwas=lapply(1:length(gwas), function(x) {
	gwas[[x]]$p_gc=pchisq(chi_gc[[x]], df=1, lower.tail=F)
	return (gwas[[x]])
	})

lapply(1:length(gwas), function(x) data.table::fwrite(
		gwas[[x]], 
		row.names=F,
		file = paste0(path_to_result_directory,'ID_',intercept_data$gwas_id[x],'_gc_corrected.csv'),
		sep = '\t')
	)

