for i in 4049 4050 4058 4059 4179
do
	run_ldscore --rg --gwas-id-1=${i} --gwas-id-2=4049,4050,4058,4059,4179 --overwrite
done

output_dir='gene_corr'
mkdir $output_dir

for i in 4049 4050 4058 4059 4179
do 
	run_ldscore_report --rg --gwas-id-1=${i} --gwas-id-2=4049,4050,4058,4059,4179 --output-dir $output_dir --output-file gene_corr_$i.txt
done
