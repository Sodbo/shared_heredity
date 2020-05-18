# Aim of this script is to estimate heritability of trats after
# shared heredity subtraction and their genetic correlations with shared heredity
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for LDL-SH
run_ldscore --h2 --gwas-id=
# 

# Estimate heritability for Triglycerides-SH
run_ldscore --h2 --gwas-id=
# 

# Estimate heritability for Cholesterol-SH
run_ldscore --h2 --gwas-id=
# 

# Estimate pairwise genetic correlations for SH and Lipids traits minus SH
run_ldscore --rg --gwas-id=,9 # sh and LDL-sh
##  (se )
run_ldscore --rg --gwas-id=,9 # sh and triglycerides-sh
##  (se )
run_ldscore --rg --gwas-id=,9 # sh and cholesterol-sh
##  (se )

# compare with the results of 08a script
