# Aim of this script is to obtain GWAS summary statistics for shared heredity

library(data.table)
library(dplyr)
source("../00_core_functions/linear_combination_v3.R") 

path_to_result_directory <- '/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/00_raw_data/'

gwas_files <- c('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/original_traits/02_unification_out/ID_1287001/ID_1287001_done.csv',
		'/mnt/polyomica/projects/shared_heredity/data/02_Lipids/original_traits/02_unification_out/ID_1287003/ID_1287003_done.csv',
		'/mnt/polyomica/projects/shared_heredity/data/02_Lipids/original_traits/02_unification_out/ID_1287004/ID_1287004_done.csv')

gwas <- lapply(gwas_files, fread)

aa <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/alphas.txt', row.names = 1)
covm <- read.table('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/pheno_corr_matrix.txt', row.names = 1, check.names = F)

rs_id <- lapply(gwas, function(x) x$rs_id)
snps <- Reduce(intersect, rs_id)
ind <- lapply(rs_id, function(x) match(snps, x))

gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])
z <- sapply(gwas_reordered, function(x) x$z)
eaf <- gwas_reordered[[1]]$eaf
sample_size <- sapply(gwas_reordered, function(x) x$n)

sh_gwas <- GWAS_linear_combination_Z_based(a = as.numeric(aa[2, ]), Z = z, covm = as.matrix(covm), N = sample_size, eaf = eaf)

sh_gwas <- mutate(sh_gwas, Z = b/se, p = pchisq(Z^2, 1, low = F))
sh_gwas <- mutate(sh_gwas, SNP = gwas_reordered[[1]]$rs_id)
sh_gwas <- mutate(sh_gwas, A1 = gwas_reordered[[1]]$ea, A2 = gwas_reordered[[1]]$ra, chr = gwas_reordered[[1]]$chr, pos = gwas_reordered[[1]]$bp,
		  			 eaf = gwas_reordered[[1]]$eaf)

head(sh_gwas, n = 2)
dir.create(path_to_result_directory)
data.table::fwrite(
	sh_gwas, 
	file = paste0(path_to_result_directory, 'SH_Lipids_GWAS.txt'),
	sep = '\t')
