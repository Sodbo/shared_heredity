# Aim of this script is to define shared hits under different thresholds

library(data.table)

setwd('/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits')

# Load extended table with clumping results

clump <- fread('./clumping_extended.txt')

# Initialize new table

threshold <- c(0.05, 1e-3, 1e-4, 1e-5, 5e-7, 5e-8) # we are going to test these thresholds
ntraits <- 4 # number of original traits + shared heredity 

shared_hits <- as.data.table(matrix(data = NA, nrow = length(threshold), ncol = ntraits + 1))
colnames(shared_hits) <- c('threshold', 'fsh', 'sh_2', 'sh_1', 'sh_0')
# fsh - number of full shared hits, SNPs significant under 5e-8 at least for one original traits and significant under defined threshold for other traits (shared heredity is not accounted here)
# sh_2 - number of hits, significant at least for one original trait under 5e-8 and under defined threshold for one more traitat
# sh_1 - number of hits, significant merely for one original trait under 5e-8 and nonsignificant for other original traits
# sh_0 - number of hits, significant under 5e-8 only for shared heredity and nonsignificant for original traits

shared_hits$threshold <- threshold

for(i in 1:length(threshold)){

	thr <- threshold[i]
	
	# Initialize columns in clumping table
	clump$fsh <- 0
	clump$sh_2 <- 0
	clump$sh_1 <- 0
	clump$sh_0 <- 0

	# select SNPs significant at 5e-8 at least for one original trait	
	not_sh0 <- which((clump$p_gc_bip <= 5e-8) | (clump$p_gc_mdd <= 5e-8) | (clump$p_gc_scz <= 5e-8))
	# SNPs that don't meet this criteria should be accounted as sh_0
	clump[-not_sh0, c('sh_0')] <- 1

	# select full shared SNPs from those which arn't sh_0
        fsh <- which((clump$p_gc_bip[not_sh0] <= thr) & (clump$p_gc_mdd[not_sh0] <= thr) & (clump$p_gc_scz[not_sh0] <= thr))
	clump$fsh[not_sh0][fsh] <- 1

	# select sh_1 SNPs
	sh1 <- which(((clump$p_gc_bip[not_sh0][-fsh] > thr) + (clump$p_gc_mdd[not_sh0][-fsh] > thr) + (clump$p_gc_scz[not_sh0][-fsh] > thr)) == 2)
	clump$sh_1[not_sh0][-fsh][sh1] <- 1

	# residual SNPs are sh_2 SNPs
	clump$sh_2[not_sh0][-fsh][-sh1] <- 1	
	
	# Fill the table
	shared_hits$sh_0[i] <- sum(clump$sh_0)
	shared_hits$sh_1[i] <- sum(clump$sh_1)
	shared_hits$sh_2[i] <- sum(clump$sh_2)
	shared_hits$fsh[i] <- sum(clump$fsh)

}	

fwrite(shared_hits, './shared_hits_different_thresholds.txt', sep = '\t')


