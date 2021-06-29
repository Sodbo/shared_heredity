# Aim of this script is to count loci significant both for the original trait and shared heredity

trait_designation <- c('bip','mdd','scz','happiness')

library(data.table)

clump <- fread('../../../data/03_neurodegenerative_diseases/several_traits/four_traits/clumping_of_SH_and_orig_traits_partIV_5e-8.txt')

ind <- grep('SH', clump$traits) # select loci significant for SH
sh <- clump[ind, ]

significant <- lapply(trait_designation, function (x) grepl(paste0('\\b', x, '\\b'), sh$traits, perl = T))
count <- sapply(significant, table)
count

