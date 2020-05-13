#Script to calculate intercept to answer whether correction is necessary.

gwas_ids=207,191,192,193,194,199 # comma separated
output_path='../../data/01_anthropometry_results/five_traits/GIP1/'
shift
echo 'Output directory: '$output_path
echo 'Gwas IDs to analyse:' $gwas_ids
run_ldscore --h2 --gwas-id=$gwas_ids --overwrite
echo 'select * from gwas.ldsc where gwas_id in ('$gwas_ids');' | psql -U gwas_user -d test_tiysgwas -h 172.25.8.65 -A -F , --pset footer -o ${output_path}intrecept_data.csv
