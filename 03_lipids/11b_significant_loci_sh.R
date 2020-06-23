# Aim of this script is to count loci significant both for the original trait and SH

trait_designation <- c(1287001, 1287003, 1287004)

library(data.table)

clump <- fread('../../../data/02_Lipids/three_traits/clumping_of_SH_and_orig_traits_partIV_5e-8.txt')

ind <- grep('SH', clump$traits) # select loci significant for SH
sh <- clump[ind, ]

significant <- lapply(trait_designation, function (x) grepl(x, sh$traits))
count <- sapply(significant, table)
colnames(count) <- trait_designation
count[2, ]

