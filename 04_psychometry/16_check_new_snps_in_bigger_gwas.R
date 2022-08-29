# Aim of this script is to check the significance of new
# (significant in SGIT but not in original traits) SNPs in a bigger GWAS on MDD

library(data.table)
library(openxlsx)

mdd <- fread("/mnt/polyomica/projects/MR_BP/00_data/RF/01_biggest_gwases/05_unified/07_Depression_done.csv", data.table = F)
snips <- read.xlsx("~/20211214_Supplementary_Tables_1-3.xlsx", sheet = "ST3b", startRow = 4)

dim(mdd)
colnames(mdd)

dim(snips)
colnames(snips)

# Find overlap
snips$p_from_biggest_gwas <- NA

i1 <- match(snips$SNP, mdd$rs_id)
table(is.na(i1))
i1 <- na.omit(i1)
mdd <- mdd[i1, ]

i2 <- match(mdd$rs_id, snips$SNP)
table(snips$SNP[i2] == mdd$rs_id)
snips$p_from_biggest_gwas[i2] <- mdd$p

fwrite(snips, "~/ST3b.txt", dec = ".", sep = "\t")

