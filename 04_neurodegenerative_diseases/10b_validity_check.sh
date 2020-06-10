# Aim of this script is to estimate heritability of trats after
# shared heredity subtraction and their genetic correlations with shared heredity
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for BIP-SH
run_ldscore --h2 --gwas-id=619
# 0.1051 ## 0.1085

# Estimate heritability for MDD-SH
run_ldscore --h2 --gwas-id=620
# 0.0591 ## 0.0576

# Estimate heritability for SCZ-SH
run_ldscore --h2 --gwas-id=621
# 0.1047 ## 0.1002

# Estimate heritability for Happiness-SH
run_ldscore --h2 --gwas-id=622
## 0.0564


# Estimate pairwise genetic correlations for SH and PGC traits minus SH
run_ldscore --rg --gwas-id=619,618 # sh and bip-sh
## -0.0071 (se 0.035)
run_ldscore --rg --gwas-id=620,618 # sh and mdd-sh
## -0.0031 (se 0.031)
run_ldscore --rg --gwas-id=621,618 # sh and scz-sh
## 0.007 (se 0.036)
run_ldscore --rg --gwas-id=622,618 # sh and happiness-sh
## 0.0055 (se 0.036)

# compare with the results of 15a script
