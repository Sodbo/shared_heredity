# Aim of this script is to obtain GWAS summary statistics for MaxH

library('data.table')
library(dplyr)
source("../00_core_functions/linear_combination_v3.R")

path_to_output_directory <- '../../../data/03_neurodegenerative_diseases/several_traits/four_traits/MaxH/00_raw_data/'
gwas_files <- c('../../../data/03_neurodegenerative_diseases/BIP/02_unification_results/bip_output_done.csv',
		'../../../data/03_neurodegenerative_diseases/MDD/02_unification_results/mdd_output_done.csv',
		'../../../data/03_neurodegenerative_diseases/SCZ/02_unification_results/scz_output_done.csv',
		'../../../data/03_neurodegenerative_diseases/happiness/02_unification_results/happiness_output_done.csv')
gwas <- lapply(gwas_files, fread)

aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/alphas.txt', row.names = 1)
covm <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)

rs_id <- lapply(gwas, function(x) x$rs_id)
snps <- Reduce(intersect, rs_id)
ind <- lapply(rs_id, function(x) match(snps, x))


gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])
z <- sapply(gwas_reordered, function(x) x$z)
sample_size <- sapply(gwas_reordered, function(x) x$n)


maxh_gwas <- GWAS_linear_combination_Z_based(a = as.numeric(aa[1, ]), Z = z, covm = as.matrix(covm), N = sample_size, eaf = gwas_reordered[[1]]$eaf)

maxh_gwas <- mutate(maxh_gwas, Z = b/se, p = pchisq(Z^2, 1, low = F))
maxh_gwas <- mutate(maxh_gwas, SNP = gwas_reordered[[1]]$rs_id)
maxh_gwas <- mutate(maxh_gwas, A1 = gwas_reordered[[1]]$ea, A2 = gwas_reordered[[1]]$ra, chr = gwas_reordered[[1]]$chr, pos = gwas_reordered[[1]]$bp,
			 eaf = gwas_reordered[[1]]$eaf)

head(maxh_gwas, n = 2)

data.table::fwrite(
	maxh_gwas, 
	row.names = F,
	file = paste0(path_to_output_directory, 'MaxH_GWAS.txt'),
	sep = '\t')
