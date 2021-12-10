# Aim of this script is to estimate heritability of UGITs
# and their genetic correlations with SGIT
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for BIP UGIT
run_ldscore --h2 --gwas-id=619
# 0.1051 ## 0.1085

# Estimate heritability for MDD UGIT
run_ldscore --h2 --gwas-id=620
# 0.0591 ## 0.0576

# Estimate heritability for SCZ UGIT
run_ldscore --h2 --gwas-id=621
# 0.1047 ## 0.1002

# Estimate heritability for Happiness UGIT
run_ldscore --h2 --gwas-id=622
## 0.0564


# Estimate pairwise genetic correlations for SGIT and PGC trait UGITs
run_ldscore --rg --gwas-id=619,618 # SGIT and bip UGIT
## -0.0071 (se 0.035)
run_ldscore --rg --gwas-id=620,618 # SGIT and mdd UGIT
## -0.0031 (se 0.031)
run_ldscore --rg --gwas-id=621,618 # SGIT and scz UGIT
## 0.007 (se 0.036)
run_ldscore --rg --gwas-id=622,618 # SGIT and happiness UGIT
## 0.0055 (se 0.036)

# compare with the results of 15a script
