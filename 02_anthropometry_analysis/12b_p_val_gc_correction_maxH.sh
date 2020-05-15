#This script makes new summary statistics file with column p_gc containing p-values corrected for genomic control

library(data.table)
path_to_maxH_directory='../../data/01_anthropometry_results/five_traits/linear_combination/'
intercept_data<-read.csv(paste0(path_to_maxH_directory,'intrecept_data_maxH.csv'))
maxH<-fread(paste0(path_to_maxH_directory,'maxH_GWAS.txt'))

chi<-qchisq(maxH$p, df=1, lower.tail=F)
chi_gc<-chi/intercept_data$intercept
maxH$p_gc<-pchisq(chi_gc, df=1, lower.tail=F)


fwrite(maxH, row.names=F,
		file = paste0(path_to_maxH_directory, 'maxH_GWAS_gc_corrected.txt'),
		sep = '\t')

