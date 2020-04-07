# Aim of this script is to estimate heritability of trats after
# shared heredity subtraction and their genetic correlations with shared heredity
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for BIP-SH
run_ldscore --h2 --gwas-id=606
# 0.1051

# Estimate heritability for MDD-SH
run_ldscore --h2 --gwas-id=609
# 0.0591

# Estimate heritability for SCZ-SH
run_ldscore --h2 --gwas-id=608
# 0.1047

# Estimate pairwise genetic correlations for SH and PGC traits minus SH
run_ldscore --rg --gwas-id=606,603 # sh and bip-sh
# -0.0066 (se 0.035)
run_ldscore --rg --gwas-id=609,603 # sh and mdd-sh
# 0.0031 (se 0.035)
run_ldscore --rg --gwas-id=608,603 # sh and scz-sh
# 0.0059 (se 0.034)

# compare with the results of 15a script
