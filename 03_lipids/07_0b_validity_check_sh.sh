# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH
run_ldscore --h2 --gwas-id=6
## 

# Estimate pairwise genetic correlations for SH and Lipids traits
run_ldscore --rg --gwas-id=3,6 # sh and choleserol
## 
run_ldscore --rg --gwas-id=4,6 # sh and LDL
## 
run_ldscore --rg --gwas-id=5,6 # sh and triglycerides
## 

# compare with the results of 07a script
