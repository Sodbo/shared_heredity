path='../../data/01_anthropometry_results/GWAS/scaled_filtered/reupload_gwases/'
sh_id=187
gwas_ids='191,192,193,194'

run_ldscore --rg --gwas-id-1=$sh_id --gwas-id-2=$(echo $gwas_ids) --overwrite

output_dir=$path'gene_corr'
mkdir $output_dir

run_ldscore_report --rg --gwas-id-1=$sh_id --gwas-id-2=$(echo $gwas_ids) --output-dir $output_dir --output-file gene_corr_$i.txt

 
