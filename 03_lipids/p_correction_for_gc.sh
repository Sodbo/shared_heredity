output_path=$1
shift
echo 'Output directory: '$output_path
echo 'Gwas IDs to analyse:' $*
ids=$(echo $* | tr ' ' ',')
run_ldscore --h2 --gwas-id=$ids --overwrite
echo 'select * from gwas.ldsc where gwas_id in ('$ids');' | psql -U gwas_user -d test_sharedheredity -h localhost -A -F , --pset footer -o ${output_path}intrecept_data.csv
