# This is an automatic pipeline to calculate correlation matrices and estimate 
# alpha coefficients for each trait in linear combination

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

# Using:
# Please, edit 00a_start.sh and run under container:
# source 00a_start.sh
{
mkdir -p $1
source 01_pheno_corr.sh $*
Rscript 02_convert_long_to_wide_form.R $1
source 03_gene_corr.sh $*
Rscript 04_gene_corr_to_matrices.R $1
Rscript 05_alpha_coefficients.R $1
Rscript 06_alpha_CI_estimation.R $1
} 2>&1 | tee $1/logs_01-06.txt
