# Aim of this script is to check the significance of new
# (significant in SGIT but not in original traits) SNPs in a bigger GWAS on BMI

library(data.table)

bmi <- fread("/mnt/polyomica/projects/shared_heredity/data/01_anthropometry_results/five_traits/GWAS/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt", data.table = F)
snips <- fread("/mnt/polyomica/projects/shared_heredity/data/01_anthropometry_results/five_traits/antro_sign_v2.txt", data.table = F)

dim(bmi)
colnames(bmi)

dim(snips)
colnames(snips)

# Prerpare rs_ids of SNPs in BMI
bmi_snp <- strsplit(bmi$"SNP", ":")
bmi_snp <- lapply(bmi_snp, function(x) x[1])
bmi_snp <- unlist(bmi_snp)

# Find overlap
snips$p_from_biggest_gwas <- NA

i1 <- match(snips$SNP, bmi_snp)
table(is.na(i1))
table(snips$SNP == bmi_snp[i1])
snips$p_from_biggest_gwas <- bmi$P[i1]

fwrite(snips, "/mnt/polyomica/projects/shared_heredity/data/01_anthropometry_results/five_traits/antro_sign_bmi_added.txt", dec = ".", sep = "\t")

