# Aim of this script is to count the loci significant both for the original trait and SH

library(data.table)

trait_designation <- c(191, 192, 193, 194, 199)

clump <- fread('../../../data/01_anthropometry_results/five_traits/joint_clumping_and_enrichment_new/clumping_of_SH_and_orig_traits_partIV_5e-8.txt')

ind <- grep('SH', clump$traits) # select loci significant for SH
sh <- clump[ind, ]

significant <- lapply(trait_designation, function (x) grepl(x, sh$traits)) # count loci significant both for SH and particular trait
count <- sapply(significant, table)
colnames(count) <- trait_designation
count[2, ]

