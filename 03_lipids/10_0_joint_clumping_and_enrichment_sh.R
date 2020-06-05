# Aim of this script is to obtain all figures and tables for SH enrichment of shared hits

library('data.table')
library('dplyr')
library('pROC')

path_for_output_results <- "/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/"

# Original traits
input_file_name <- c('/mnt/polyomica/projects/shared_heredity/data/02_Lipids/original_traits/02_unification_out/ID_1287001/ID_1287001_done.csv',
		     '/mnt/polyomica/projects/shared_heredity/data/02_Lipids/original_traits/02_unification_out/ID_1287003/ID_1287003_done.csv',
		     '/mnt/polyomica/projects/shared_heredity/data/02_Lipids/original_traits/02_unification_out/ID_1287004/ID_1287004_done.csv')
gwas <- lapply(input_file_name, fread)
gwas <- lapply(gwas, function(x) x[(x$p != 0) & (x$eaf >= 0.01) & (x$eaf <= 0.99), ])

# GWAS for shared heredity
gwas_sh <- fread("/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/02_unification_results/sh_output_done.csv", data.table = F)
gwas_sh <- gwas_sh[(gwas_sh$p != 0) & (gwas_sh$eaf >= 0.01) & (gwas_sh$eaf <= 0.99), ]

# Reordering of all original GWASs and SH
rs_id <- lapply(gwas, function(x) x$rs_id)
snps <- Reduce(intersect, rs_id)
snps <- intersect(snps, gwas_sh$rs_id)

ind <- lapply(rs_id, function(x) match(snps, x))
gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])

ind <- match(snps, gwas_sh$rs_id)
gwas_sh <- gwas_sh[ind, ]

# Forming the input tables 
pval_orig <- sapply(gwas_reordered, function(x) x$p)
pval_sh <- gwas_sh$p
chr <- gwas_sh$chr
pos <- gwas_sh$bp


source("../00_core_functions/joint_function_for_enrichment_and_auc.R")
	
out <- joint_function(pval_orig = pval_orig, pval_sh = pval_sh, chr = chr, pos = pos, thr_sh_hits = 1e-5, thr_sh_sig = 5e-8,
	orig_traits_names = input_file_name, path_output = path_for_output_results, delta = 5e5, SNPs = snps)
