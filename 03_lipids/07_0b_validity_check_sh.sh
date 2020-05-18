# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH
run_ldscore --h2 --gwas-id=9
## 

# Estimate pairwise genetic correlations for SH and Lipids traits
run_ldscore --rg --gwas-id=6,9 # sh and LDL
## 
run_ldscore --rg --gwas-id=7,9 # sh and triglycerides
## 
run_ldscore --rg --gwas-id=8,9 # sh and choleserol
## 

# compare with the results of 07a script
