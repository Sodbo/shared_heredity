# Aim of this script is to make correction for Genomic control (GC) for SGIT

library(data.table)

setwd('/mnt/polyomica/projects/shared_heredity/')

gwas_file <- './data/02_Lipids/three_traits/SH/02_unification_results/sh_output_done.csv.gz'

intercept_data <- read.csv('./data/02_Lipids/traits_minus_SH/three_traits/gene_corr/h2.csv') # load table with parameters for GC 
gwas <- fread(gwas_file)

# Make GC correction for intercept
chi <- qchisq(gwas$p, df = 1, lower.tail = F)
chi_gc <- chi/intercept_data$intercept[intercept_data$gwas_id == 9]
gwas$p_gc <- pchisq(chi_gc, df = 1, lower.tail = F)

# Write corrected GWAS summary statistics
fwrite(gwas, file = './data/02_Lipids/three_traits/SH/03_gc_corrected/sh_gc_corrected.csv', sep = '\t')


