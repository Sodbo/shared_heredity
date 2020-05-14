# Aim of this script is to make correction for Genomic control (GC) for original traits and shared heredity

library(data.table)

intercept_data <- read.csv('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/GIP/intrecept_data.csv') # load table with parameters for GC 
gwas <- fread('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/GIP/02_unification_results/gip1_output_done.csv')

# Make GC correction for intercept
chi <- qchisq(gwas$p, df = 1, lower.tail = F)
chi_gc <- chi/intercept_data$intercept[1]
gwas$p_gc <- pchisq(chi_gc, df = 1, lower.tail = F)

# Write corrected GWAS summary statistics
fwrite(gwas, row.names=F, file = '../../../data/03_neurodegenerative_diseases/several_traits/four_traits/GIP/03_gc_corrected/gip1_gc_corrected.csv', sep = '\t')

