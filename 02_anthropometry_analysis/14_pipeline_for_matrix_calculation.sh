# Path to output directory should be given as the first command line argument
# Ids of GWAS should also be given as the command line arguments after path.
# Necessary R libraries should be installed before the first run.
# Results will be written to the following files:
# gene_corr_matrix.txt
# gene_corr directory
# gene_cov_matrix.txt
# phen_corr_res.txt
# pheno_corr_matrix.txt
# alphas.txt
# w.txt

# Example of using:
# source start.sh

mkdir $1
source 01_pheno_corr.sh $*
Rscript 02_convert_long_to_wide_form.R $1
source 03_gene_corr.sh $*
Rscript 04_gene_corr_to_matrix.R $1
Rscript 05_shared_heredity.R $1
