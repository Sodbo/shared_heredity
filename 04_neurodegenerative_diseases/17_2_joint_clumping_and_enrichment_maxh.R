# Aim of this script is to obtain all figures and tables for MaxH enrichment of shared hits

library('data.table')
library('dplyr')
library('pROC')

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases")

path_for_output_results <- "../../../data/03_neurodegenerative_diseases/several_traits/four_traits/MaxH/"

#GC corrected original GWASes
input_file_name <- c('../../../data/03_neurodegenerative_diseases/BIP/03_gc_corrected/bip_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/MDD/03_gc_corrected/mdd_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/SCZ/03_gc_corrected/scz_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/happiness/03_gc_corrected/happiness_gc_corrected.csv')
gwas <- lapply(input_file_name, fread)

#GWAS for MaxH corrected for GC
gwas_maxh <- fread("../../../data/03_neurodegenerative_diseases/several_traits/four_traits/MaxH/03_gc_corrected/maxh_gc_corrected.csv", data.table = F)

#reordering of all original GWASs and MaxH
rs_id <- lapply(gwas, function(x) x$rs_id)
snps <- Reduce(intersect, rs_id)
snps <- intersect(snps, gwas_maxh$rs_id)

ind <- lapply(rs_id, function(x) match(snps, x))
gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])

ind <- match(snps, gwas_maxh$rs_id)
gwas_maxh <- gwas_maxh[ind, ]

# Forming the input tables 
pval_orig <- sapply(gwas_reordered, function(x) x$p_gc)
pval_maxh <- gwas_maxh$p_gc
chr <- gwas_maxh$chr
pos <- gwas_maxh$bp


source("../00_core_functions/joint_function_for_enrichment_and_auc.R")
	
out <- joint_function(pval_orig = pval_orig, pval_sh = pval_maxh, chr = chr, pos = pos, thr_sh_hits = 1e-3, thr_sh_sig = 5e-8,
	orig_traits_names = input_file_name, path_output = path_for_output_results, delta = 5e5, SNPs = snps)
