# Aim of this script is to perform joint clumping for original traits, SGCT and UGCTs

library('data.table')
library('dplyr')

setwd("/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_psychometry/")

path_for_output_results <- "../../../data/03_neurodegenerative_diseases/several_traits/four_traits/joint_clumping_full_set/"
dir.create(path_for_output_results)

#GC corrected GWASes on original traits (except MDD - it is not corrected) and UGCTs
input_file_name <- c('../../../data/03_neurodegenerative_diseases/BIP/03_gc_corrected/bip_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/MDD/02_unification_results/mdd_output_done.csv.gz',
		     '../../../data/03_neurodegenerative_diseases/SCZ/03_gc_corrected/scz_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/happiness/03_gc_corrected/happiness_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/BIP/03_gc_corrected/bip-sh_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/MDD/03_gc_corrected/mdd-sh_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/SCZ/03_gc_corrected/scz-sh_gc_corrected.csv',
		     '../../../data/03_neurodegenerative_diseases/happiness/03_gc_corrected/happiness-sh_gc_corrected.csv')
gwas <- lapply(input_file_name, fread)

#GWAS for shared heredity corrected for GC
gwas_sh <- fread("../../../data/03_neurodegenerative_diseases/several_traits/four_traits/linear_combination/03_gc_corrected/sh_gc_corrected.csv", data.table = F)

#reordering of all original GWASs and SH
rs_id <- lapply(gwas, function(x) x$rs_id)
snps <- Reduce(intersect, rs_id)
snps <- intersect(snps, gwas_sh$rs_id)

ind <- lapply(rs_id, function(x) match(snps, x))
gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])

ind <- match(snps, gwas_sh$rs_id)
gwas_sh <- gwas_sh[ind, ]

# Forming the input tables
pval_bip <- as.matrix(gwas_reordered[[1]]$p_gc)
pval_mdd <- as.matrix(gwas_reordered[[2]]$p)
pval_orig <- sapply(gwas_reordered[3:8], function(x) x$p_gc)
pval_orig <- cbind(pval_bip, pval_mdd, pval_orig)
pval_sh <- gwas_sh$p_gc
chr <- gwas_sh$chr
pos <- gwas_sh$bp

trait_names<-c('BIP', 'MDD', 'SCZ', 'Happiness', 'BIP UGCT', 'MDD UGCT', 'SCZ UGCT', 'Happiness UGCT')

source("../00_core_functions/joint_function_for_enrichment_and_auc.R")
	
out <- joint_function(pval_orig = pval_orig, pval_sh = pval_sh, chr = chr, pos = pos, thr_sh_hits = 1e-3, thr_sh_sig = 5e-8,
	orig_traits_names = trait_names, path_output = path_for_output_results, delta = 5e5, SNPs = snps)
