# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH
run_ldscore --h2 --gwas-id=618
## 0.4075 (for 4 traits)

# Estimate pairwise genetic correlations for SH and PGC traits
run_ldscore --rg --gwas-id=551,618 # sh and bip
## 0.9415 (for 4 traits)
run_ldscore --rg --gwas-id=554,618 # sh and mdd
## 0.4926 (for 4 traits)
run_ldscore --rg --gwas-id=555,618 # sh and scz
## 0.8964 (for 4 traits)
run_ldscore --rg --gwas-id=617,618 # sh and happiness
## 0.359

# compare with the results of 07a script
