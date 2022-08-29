# Aim of this script is to check the significance of new
# (significant in SGIT but not in original traits) SNPs in a bigger GWAS on BMI

library(data.table)
library(openxlsx)

bmi <- fread("/mnt/polyomica/projects/shared_heredity/data/01_anthropometry_results/five_traits/GWAS/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt", data.table = F)
#mdd <- fread("/mnt/polyomica/projects/MR_BP/00_data/RF/01_biggest_gwases/05_unified/07_Depression_done.csv, data.table = F)
#snips <- read.xlsx("~/20211214_Supplementary_Tables_1-3.xlsx", sheet = "ST3b", startRow = 4)
snips <- read.xlsx("~/20211214_Supplementary_Tables_1-3.xlsx", sheet = "ST3a", startRow = 4)

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

fwrite(snips, "~/ST3a.txt", dec = ".", sep = "\t")

