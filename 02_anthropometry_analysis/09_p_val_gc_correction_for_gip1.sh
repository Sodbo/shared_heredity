#This script makes new summary statistics file with column p_gc containing p-values corrected for genomic control

library(data.table)
path_to_gip1_directory='../../data/01_anthropometry_results/five_traits/GIP1/'
intercept_data<-read.csv(paste0(path_to_gip1_directory,'intrecept_data.csv'))
gip1<-fread(paste0(path_to_gip1_directory,'GIP1_GWAS.txt'))

chi<-qchisq(gip1$p, df=1, lower.tail=F)
chi_gc<-chi/intercept_data$intercept
gip1$p_gc<-pchisq(chi_gc, df=1, lower.tail=F)


fwrite(gip1, row.names=F,
		file = paste0(path_to_gip1_directory, 'GIP1_GWAS_gc_corrected.txt'),
		sep = '\t')

