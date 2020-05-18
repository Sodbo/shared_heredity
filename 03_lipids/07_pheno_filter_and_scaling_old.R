# Aim of this script is to obtain GWAS summary statistics with variance scaled to 1

library(dplyr)
library(data.table)
sd.table <-fread(
	file = '/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/sd.txt',
	header = TRUE,
	data.table=F
)

gwas.ids=colnames(sd.table)
dir.create('../data/02_Lipids/GWAS/scaled_filtered/', showWarnings = FALSE)
for (col.number in 1:ncol(sd.table)) {
	sd.trait <- sd.table[, col.number]
	trait<-gwas.ids[col.number]
	cat(trait,"; ","SD: ",sd.trait,"\n")	
	# Reading raw GWAS-file

	gwas.name <- paste0('ID_', trait, '.csv')
	raw.gwas <- fread(
		input = paste0('/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/GWAS/',gwas.name),
		data.table=F,
		  header=T,
		  stringsAsFactors=F)

	#gwas.info.filtered <- filter(raw.gwas,info >= 0.7) # Filtering by info
	#please, change 'raw.gwas' with 'gwas.info.filtered' in the next expression if you uncomment the above expression
	
	gwas.info.filtered=mutate(raw.gwas, MAF=pmin(eaf,1-eaf))
	gwas.MAF.filtered <- filter(gwas.info.filtered, MAF >= 1e-5) # Filtering by info
	
	cat("Nsnps after filtering:", nrow(gwas.MAF.filtered),"\n")
	gwas.standart <- gwas.MAF.filtered
	gwas.standart$beta <- gwas.standart$beta / sd.trait
	gwas.standart$se <- gwas.standart$se / sd.trait
	
	# Writing an output file
	fwrite(gwas.standart, 
		file = paste0('/home/ubuntu/polyomica/projects/shared_heredity/data/02_Lipids/GWAS/scaled_filtered/', gwas.name))

}


