# Aim of this script is to estimate heritability of
# GIP1 and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for GIP1
run_ldscore --h2 --gwas-id=623
## 0.4126

# Estimate pairwise genetic correlations for GIP1 and PGC traits
run_ldscore --rg --gwas-id=551,623 # gip1 and bip
## 0.9485
run_ldscore --rg --gwas-id=554,623 # gip1 and mdd
## 0.4724
run_ldscore --rg --gwas-id=555,623 # gip1 and scz
## 0.8931
run_ldscore --rg --gwas-id=617,623 # gip1 and happiness
## 0.3296

# compare with the results of 07_1a script
