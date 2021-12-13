# Aim of this script is to estimate heritability of UGITs
# and their genetic correlations with SGIT
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for LDL UGIT
run_ldscore --h2 --gwas-id=10
# 0.0994 

# Estimate heritability for Triglycerides UGIT
run_ldscore --h2 --gwas-id=11
# 0.1714

# Estimate heritability for Cholesterol UGIT
run_ldscore --h2 --gwas-id=12
# 0.1079

# Estimate pairwise genetic correlations for SGIT and Lipid trait UGITs
run_ldscore --rg --gwas-id=10,9 # SGIT and LDL UGIT
## -0.4904 (se 0.693)
run_ldscore --rg --gwas-id=11,9 # SGIT and triglycerides UGIT
## 0.1671 (se 0.27)
run_ldscore --rg --gwas-id=12,9 # SGIT and cholesterol UGIT
## 0.0805 (se 0.106)

# compare with the results of 08a script
