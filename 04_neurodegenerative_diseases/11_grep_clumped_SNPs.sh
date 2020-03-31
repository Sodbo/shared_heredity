# This script is to grep clumped SNPs from original traits and shared heredity

path='/mnt/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases'
while IFS= read -r line; do
	grep -w "$line" $path'/BIP/03_gc_corrected/bip_gc_corrected.csv' >> $path'/BIP/04_clumped_snps/bip_clumped_snps_only.txt'
	grep -w "$line" $path'/MDD/03_gc_corrected/mdd_gc_corrected.csv' >> $path'/MDD/04_clumped_snps/mdd_clumped_snps_only.txt'
	grep -w "$line" $path'/SCZ/03_gc_corrected/scz_gc_corrected.csv' >> $path'/SCZ/04_clumped_snps/scz_clumped_snps_only.txt'
	grep -w "$line" $path'/several_traits/linear_combination/03_gc_corrected/sh_gc_corrected.csv' >> $path'/several_traits/linear_combination/04_clumped_snps/sh_clumped_snps_only.txt'
done < $path'/several_traits/Clumped_SNPs_original_traits_and_SH.txt'
