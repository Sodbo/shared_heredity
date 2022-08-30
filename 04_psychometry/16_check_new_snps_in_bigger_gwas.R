# Aim of this script is to check the significance of new
# (significant in SGIT but not in original traits) SNPs in a bigger GWAS on MDD

library(data.table)

mdd <- fread("/mnt/polyomica/projects/MR_BP/00_data/RF/01_biggest_gwases/05_unified/07_Depression_done.csv", data.table = F)
snips <- fread("/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/four_traits/psy_sign_v2.txt", data.table = F)

dim(mdd)
colnames(mdd)

dim(snips)
colnames(snips)

# Find overlap
snips$p_from_biggest_gwas <- NA

i1 <- match(snips$SNP, mdd$rs_id)
which(is.na(i1))

snips$p_from_biggest_gwas <- mdd$p[i1]
table(snips$p_from_biggest_gwas == mdd$p[i1])

fwrite(snips, "/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/four_traits/psy_sign_mdd_added.txt", dec = ".", sep = "\t")

