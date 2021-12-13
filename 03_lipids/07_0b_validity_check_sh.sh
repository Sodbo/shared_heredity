# Aim of this script is to estimate heritability of
# SGIT and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH
run_ldscore --h2 --gwas-id=9
## h2=0.2126  se=0.0295

# Estimate pairwise genetic correlations for SGIT and Lipids traits
run_ldscore --rg --gwas-id=6,9 # SGIT and LDL
## 0.9605 se=0.02
run_ldscore --rg --gwas-id=7,9 # SGIT and triglycerides
## 0.7404 se=0.129
run_ldscore --rg --gwas-id=8,9 # SGIT and choleserol
## 1.0177 se=0.018

# compare with the results of 07a script
