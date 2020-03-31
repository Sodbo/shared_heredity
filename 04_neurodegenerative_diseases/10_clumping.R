# Aim of this script is to do clumping for original traits and SH and generate a final table

library(data.table)

setwd('/mnt/polyomica/projects/shared_heredity/elgaeva_src/shared_heredity/04_neurodegenerative_diseases/')
source("../00_core_functions/clumping_v2.R", chdir = F)

# Define files with summary statistics of original traits and shared heredity after genomic control
input_file_name<-c('../../../data/03_neurodegenerative_diseases/BIP/03_gc_corrected/bip_gc_corrected.csv',
		   '../../../data/03_neurodegenerative_diseases/MDD/03_gc_corrected/mdd_gc_corrected.csv',
		   '../../../data/03_neurodegenerative_diseases/SCZ/03_gc_corrected/scz_gc_corrected.csv',
		   '../../../data/03_neurodegenerative_diseases/several_traits/linear_combination/03_gc_corrected/sh_gc_corrected.csv')

thr <- 5e-8 # set variable

out <- NULL # initialize variable

# Generate a table with clumping results for original traits and shared heredity

for (input in input_file_name) {

	sst <- fread(input, header = T, stringsAsFactors = F, data.table = F)
	sst_sm <- sst[(pmin(sst$eaf, 1 - sst$eaf) >= 0.01), ] # maf filter
	trait <- strsplit(input, '03_gc_corrected/')[[1]][2]
	trait <- strsplit(trait, '_')[[1]][1]
	
	lt <- function_for_shlop_29_03_2020(sst_sm, p_value = "p_gc", pos = "bp", snp = "rs_id", delta = 5e5, chr = "chr", thr = thr)
	
	if (nrow(lt) > 0 ) {
		lt <- cbind(lt, trait)
		out <- rbind(out, lt)
	}
}

dim(out)
table(out$trait)
#bip mdd scz  sh
# 14   3  97  48


# Repeat clumping procedure and count signals from the same locus for different traits
final <- function_for_shlop_29_03_2020(out, trait = "trait", p_value = "p_gc", pos = "bp", snp = "rs_id", delta = 5e5, chr = "chr")
final <- final[order(final$chr, final$bp), ]
colnames(final)
str(final)
final$trait <- as.character(final$trait)
str(final)

dim(final)
colnames(final)

# Write table with clumping results for original traits and shared heredity
path <- '../../../data/03_neurodegenerative_diseases/several_traits/'
result_file_name <- 'Clumping_for_all_orignal_traits_and_SH.txt'

fwrite(final, paste0(path, result_file_name), row.names = F, sep = '\t', dec = '.')

# Write file with unique clumped SNPs for original traits and shared heredity
write.table(final$rs_id, file = paste0(path, 'Clumped_SNPs_original_traits_and_SH.txt'), row.names = F, col.names = F, quote = F)

