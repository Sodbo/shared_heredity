library(data.table)
path_to_result_directory<-'../../data/anthropometry_results/four_traits/GWAS/scaled_filtered/'
gwas_files<-list.files(path_to_result_directory, full.names=T, pattern="ID_\\d+.csv")
intercept_data<-read.csv('../../data/anthropometry_results/four_traits/intrecept_data.csv') 
gwas<-lapply(gwas_files, fread)
chi=lapply(gwas, function(x) qchisq(x$p, df=1, lower.tail=F))

chi_gc=lapply(1:length(chi), function(x) chi[[x]]/intercept_data$intercept[x])
gwas=lapply(1:length(gwas), function(x) {
	gwas[[x]]$p_gc=pchisq(chi_gc[[x]], df=1, lower.tail=F)
	return (gwas[[x]])
	})

lapply(1:length(gwas), function(x) data.table::fwrite(
		gwas[[x]], 
		row.names=F,
		file = paste0('../../data/anthropometry_results/four_traits/GWAS/scaled_filtered/ID_',intercept_data$gwas_id[x],'_gc_corrected.csv'),
		sep = '\t')
	)

