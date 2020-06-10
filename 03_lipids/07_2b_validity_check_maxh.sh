# Aim of this script is to estimate heritability of
# MaxH and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for MaxH
run_ldscore --h2 --gwas-id=14
## 0.214

# Estimate pairwise genetic correlations for MaxH and Lipid traits
run_ldscore --rg --gwas-id=6,14 # maxh and LDL
## 0.7766
run_ldscore --rg --gwas-id=7,14 # maxh and triglycerides
## 0.9603
run_ldscore --rg --gwas-id=8,14 # maxh and cholesterol
## 0.8074

