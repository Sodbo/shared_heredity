# Aim of this script is to make correction for Genomic control (GC) for original traits and shared heredity

library(data.table)

gwas_files <- c('../../../data/03_neurodegenerative_diseases/BIP/02_unification_results/bip_output_done.csv',
	      '../../../data/03_neurodegenerative_diseases/MDD/02_unification_results/mdd_output_done.csv',
	      '../../../data/03_neurodegenerative_diseases/SCZ/02_unification_results/scz_output_done.csv',
	      '../../../data/03_neurodegenerative_diseases/several_traits/linear_combination/02_unification_results/sh_output_done.csv')

intercept_data <- read.csv('../../../data/03_neurodegenerative_diseases/several_traits/intrecept_data.csv') # load table with parameters for GC 
gwas <- lapply(gwas_files, fread)

# Make GC correction for intercept
chi <- lapply(gwas, function(x) qchisq(x$p, df = 1, lower.tail = F))
chi_gc <- lapply(1:length(chi), function(x) chi[[x]]/intercept_data$intercept[x])
gwas <- lapply(1:length(gwas), function(x) {
	gwas[[x]]$p_gc <- pchisq(chi_gc[[x]], df = 1, lower.tail = F)
	return (gwas[[x]])
	})

# Define trait names
traits <- c('bip', 'mdd', 'scz', 'sh')

# Write corrected GWAS summary statistics
lapply(1:length(gwas), function(x) data.table::fwrite(
		gwas[[x]], 
		row.names=F,
		file = paste0(strsplit(gwas_files[x], '02')[[1]][1], '/03_gc_corrected/', traits[x], '_gc_corrected.csv'),
		sep = '\t')
	)

