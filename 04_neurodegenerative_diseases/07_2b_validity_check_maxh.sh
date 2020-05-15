# Aim of this script is to estimate heritability of
# MaxH and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for MaxH
run_ldscore --h2 --gwas-id=624
## 0.4165

# Estimate pairwise genetic correlations for MaxH and PGC traits
run_ldscore --rg --gwas-id=551,624 # maxh and bip
## 0.9565
run_ldscore --rg --gwas-id=554,624 # maxh and mdd
## 0.4541
run_ldscore --rg --gwas-id=555,624 # maxh and scz
## 0.8842
run_ldscore --rg --gwas-id=617,624 # maxh and happiness
## 0.3145

# compare with the results of 07_2a script
