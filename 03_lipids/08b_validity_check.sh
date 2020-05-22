# Aim of this script is to estimate heritability of trats after
# shared heredity subtraction and their genetic correlations with shared heredity
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for LDL-SH
run_ldscore --h2 --gwas-id=10
# 0.0994 

# Estimate heritability for Triglycerides-SH
run_ldscore --h2 --gwas-id=11
# 0.1714

# Estimate heritability for Cholesterol-SH
run_ldscore --h2 --gwas-id=12
# 0.1079

# Estimate pairwise genetic correlations for SH and Lipids traits minus SH
run_ldscore --rg --gwas-id=10,9 # sh and LDL-sh
## -0.4904 (se 0.693)
run_ldscore --rg --gwas-id=11,9 # sh and triglycerides-sh
## 0.1671 (se 0.27)
run_ldscore --rg --gwas-id=12,9 # sh and cholesterol-sh
## 0.0805 (se 0.106)

# compare with the results of 08a script
