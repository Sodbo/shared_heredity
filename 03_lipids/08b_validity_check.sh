# Aim of this script is to estimate heritability of trats after
# shared heredity subtraction and their genetic correlations with shared heredity
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for LDL-SH
run_ldscore --h2 --gwas-id=10
# h2=0.0994  (se=0.0198)

# Estimate heritability for Triglycerides-SH
run_ldscore --h2 --gwas-id=11
# h2=0.1714  (se=0.0257)

# Estimate heritability for Cholesterol-SH
run_ldscore --h2 --gwas-id=12
# h2=0.1079  (se=0.0185)

# Estimate pairwise genetic correlations for SH and Lipids traits minus SH
run_ldscore --rg --gwas-id=10,9 # sh and LDL-sh
## -0.4904  (se=0.68)
run_ldscore --rg --gwas-id=11,9 # sh and triglycerides-sh
## 0.1671  (se=0.28)
run_ldscore --rg --gwas-id=12,9 # sh and cholesterol-sh
## 0.0805  (se=0.106)

# compare with the results of 08a script
