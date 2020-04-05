# This script is to grep clumped SNPs from original traits and shared heredity

path='../../data/01_anthropometry_results/GWAS/scaled_filtered/shared_hits/'
grep -w -f $path'list_of_clumped_rs_id.txt' $path'../ID_4049_gc_corrected.csv' >> $path'ID_4049_clumped_snps_v2.txt'
grep -w -f $path'list_of_clumped_rs_id.txt' $path'../ID_4050_gc_corrected.csv' >> $path'ID_4050_clumped_snps_v2.txt'
grep -w -f $path'list_of_clumped_rs_id.txt' $path'../ID_4058_gc_corrected.csv' >> $path'ID_4058_clumped_snps_v2.txt'
grep -w -f $path'list_of_clumped_rs_id.txt' $path'../ID_4179_gc_corrected.csv' >> $path'ID_4179_clumped_snps_v2.txt'
grep -w -f $path'list_of_clumped_rs_id.txt' '../../data/01_anthropometry_results/linear_combination/02_unification_out/SH/SH_done.csv' >> $path'SH_clumped_snps_v2.txt'

