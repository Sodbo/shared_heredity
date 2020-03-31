# Aim of this script is to make Fisher test for Shared Heredity

library(data.table)

setwd('/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits')

# Load extended table with clumping results

clump <- fread('./clumping_extended.txt', stringsAsFactors = F) 

# Initialize new columns

clump$all_significant_and_SH <- 0
clump$not_significant_and_SH <- 0 
clump$all_significant_not_SH <- 0
clump$not_significant_not_SH <- 0 

clump$is_SH <- grepl('\\bsh\\b', clump$traits, perl = T) # SH is significant

# Redefine new columns

clump$all_significant_and_SH[clump$Ntraits == 4] <- 1 # all traits are significant including SH
clump$not_significant_and_SH[clump$Ntraits != 4 & clump$is_SH == T] <- 1 # not all traits are significant but SH is significant
clump$all_significant_not_SH[clump$Ntraits == 4 & clump$is_SH == F] <- 1 # all traits are significant except SH
clump$not_significant_not_SH[clump$Ntraits != 4 & clump$is_SH == F] <- 1 # not all significant and SH is not significant too

# Count values for table

n11 <- sum(clump$all_significant_and_SH)
n12 <- sum(clump$not_significant_and_SH)
n21 <- sum(clump$all_significant_not_SH)
n22 <- sum(clump$not_significant_not_SH)

contigency_table <- matrix(c(n11, n12, n21, n22), nrow = 2, ncol = 2, byrow = T)

colnames(contigency_table) <- c('all_significant', 'not_all_significant')
rownames(contigency_table) <- c('SH_is_significant', 'SH_is_not_significant')

write.table(contigency_table, './contigency_table.txt')

result <- fisher.test(contigency_table)

write.table(result$p.value, './SH_enrichment_analysis.txt')
