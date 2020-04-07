# This script is to substract shared genetic component from the original traits
# and to check for small correlation between traits without shared genetic component

library(data.table)
library(dplyr)

setwd('/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases')

source("../00_core_functions/linear_combination_v3.R")
source("../00_core_functions/gcov_for_linear_comb_with_i_trait.R")
source("../00_core_functions/heritability_of_linear_combination.R")

gwas_files <- c('../../../data/03_neurodegenerative_diseases/BIP/02_unification_results/bip_output_done.csv',
		'../../../data/03_neurodegenerative_diseases/MDD/02_unification_results/mdd_output_done.csv',
		'../../../data/03_neurodegenerative_diseases/SCZ/02_unification_results/scz_output_done.csv')

gwas <- lapply(gwas_files, fread)

aa <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/alphas.txt', row.names = 1)
phem <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/pheno_corr_matrix.txt', row.names = 1,  check.names = F)
gcov <- read.table('../../../data/03_neurodegenerative_diseases/several_traits/gene_cov_matrix.txt', row.names = 1,  check.names = F)

rs_id <- lapply(gwas, function(x) x$rs_id)
snps <- Reduce(intersect, rs_id)
ind <- lapply(rs_id, function(x) match(snps, x))

gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])
z <- sapply(gwas_reordered, function(x) x$z)
eaf <- gwas_reordered[[1]]$eaf
sample_size <- sapply(gwas_reordered, function(x) x$n)


alphas <- as.numeric(aa[2, ])
n_traits <- c(1:length(alphas))
slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))
position <- diag(length(alphas))


tr_sh_gwas <- lapply(n_traits, function(x) GWAS_linear_combination_Z_based(a = position[x, ] - alphas*slope[x], Z = z, covm = as.matrix(phem), N = sample_size, eaf = eaf))
tr_sh_gwas <- lapply(tr_sh_gwas, function(x) mutate(x, Z = b/se, p = pchisq(Z^2, 1, low = F)))
tr_sh_gwas <- lapply(tr_sh_gwas, function(x) mutate(x, SNP = gwas_reordered[[1]]$rs_id))
tr_sh_gwas <- lapply(tr_sh_gwas, function(x) mutate(x, A1 = gwas_reordered[[1]]$ea, A2 = gwas_reordered[[1]]$ra, chr = gwas_reordered[[1]]$chr, pos = gwas_reordered[[1]]$bp, eaf = gwas_reordered[[1]]$eaf))

lapply(tr_sh_gwas, function(x) head(x, n = 2))

# Define trait names
traits <- c('bip', 'mdd', 'scz')

lapply(n_traits, function(x) data.table::fwrite(
		tr_sh_gwas[[x]], 
		row.names = F,
		file = paste0(strsplit(gwas_files[x], '02')[[1]][1], '/05_minus_sh/', traits[x],'-sh_gwas.txt'),
		sep = '\t')
	)

