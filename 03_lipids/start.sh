# Ids of GWAS should be given as command line arguments.
# Necessary R libraries should be installed before the first run.
# Results will be written to the following files:
# gene_corr_matrix.txt
# gene_corr directory
# gene_cov_matrix.txt
# phen_corr_res.txt
# pheno_corr_matrix.txt
# alphas.txt
# w.txt
 
source 01_pheno_corr.sh $*
Rscript 02_convert_long_to_wide_form.R
source 03_gene_corr.sh $*
Rscript 04_gene_corr_to_matrix.R
#Rscript 05_shared_heredity.R
