#The script for grep clumped original traits SNPs and shared heredity SNPs

path='../../data/01_anthropometry_results/GWAS/scaled_filtered/'
for rs in $(cat ${path}list_of_clumped_rs_id.txt) 
do
echo Searching for rs_id=$rs
	grep -P "$rs\t" $path'ID_4049_gc_corrected.csv' >> $path'ID_4049_clumped_with_SH.txt'
	grep -P "$rs\t" $path'ID_4050_gc_corrected.csv' >> $path'ID_4050_clumped_with_SH.txt'
	grep -P "$rs\t" $path'ID_4058_gc_corrected.csv' >> $path'ID_4058_clumped_with_SH.txt'
	grep -P "$rs\t" $path'ID_4179_gc_corrected.csv' >> $path'ID_4179_clumped_with_SH.txt'
	grep -P "$rs\t" '../../data/01_anthropometry_results/linear_combination/SH_GWAS.txt' >> $path'SH_clumped_with_traits.txt'
done
