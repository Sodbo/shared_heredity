# Aim of this script is to estimate heritability of
# GIP1 and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for GIP1
run_ldscore --h2 --gwas-id=13
## 0.2143

# Estimate pairwise genetic correlations for GIP1 and Lipid traits
run_ldscore --rg --gwas-id=6,13 # gip1 and LDL
## 0.9154
run_ldscore --rg --gwas-id=7,13 # gip1 and triglycerides
## 0.8247
run_ldscore --rg --gwas-id=8,13 # gip1 and cholesterol
## 0.9655


# compare with the results of 07_1a script
