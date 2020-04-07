# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH
run_ldscore --h2 --gwas-id=603
# 0.4135

# Estimate pairwise genetic correlations for SH and PGC traits
run_ldscore --rg --gwas-id=551,603 # sh and bip
# 0.9521
run_ldscore --rg --gwas-id=554,603 # sh and mdd
# 0.4523
run_ldscore --rg --gwas-id=555,603 # sh and scz
# 0.8938

# compare with the results of 07a script
