# Aim of this script is to make correction for Genomic control (GC) for original traits, SGIT and UGITs

library(data.table)

gwas_files <- c('../../../data/03_neurodegenerative_diseases/BIP/02_unification_results/bip_output_done.csv',
	      '../../../data/03_neurodegenerative_diseases/MDD/02_unification_results/mdd_output_done.csv',
	      '../../../data/03_neurodegenerative_diseases/SCZ/02_unification_results/scz_output_done.csv',
	      '../../../data/03_neurodegenerative_diseases/happiness/02_unification_results/happiness_output_done.csv',
	      '../../../data/03_neurodegenerative_diseases/several_traits/four_traits/linear_combination/02_unification_results/sh_output_done.csv',
	      '/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/BIP/07_minus_sh_unified/bip-sh_output_done.csv',
	      '/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/MDD/07_minus_sh_unified/mdd-sh_output_done.csv',
	      '/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/SCZ/07_minus_sh_unified/scz-sh_output_done.csv',
	      '/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/07_minus_sh_unified/happiness-sh_output_done.csv')

intercept_data <- read.csv('/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/traits_minus_SH/four_traits/gene_corr/h2.csv') # load table with parameters for GC 
gwas <- lapply(gwas_files, fread)

# Make GC correction for intercept
chi <- lapply(gwas, function(x) qchisq(x$p, df = 1, lower.tail = F))
chi_gc <- lapply(1:length(chi), function(x) chi[[x]]/intercept_data$intercept[x])
gwas <- lapply(1:length(gwas), function(x) {
	gwas[[x]]$p_gc <- pchisq(chi_gc[[x]], df = 1, lower.tail = F)
	return (gwas[[x]])
	})

# Define trait names
traits <- c('bip', 'mdd', 'scz', 'happiness', 'sh', 'bip-sh', 'mdd-sh', 'scz-sh', 'happiness-sh')

# Write corrected GWAS summary statistics
lapply(1:length(gwas), function(x) 
       		if(grepl('02', gwas_files[x]))
       		data.table::fwrite(
		gwas[[x]], 
		row.names=F,
		file = paste0(strsplit(gwas_files[x], '02')[[1]][1], '/03_gc_corrected/', traits[x], '_gc_corrected.csv'),
		sep = '\t')
	else
		data.table::fwrite(
		gwas[[x]],
		row.names=F,
		file = paste0(strsplit(gwas_files[x], '07')[[1]][1], '/03_gc_corrected/', traits[x], '_gc_corrected.csv'),
                sep = '\t')
	)

