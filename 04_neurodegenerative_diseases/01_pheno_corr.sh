output_path=$1
shift
echo 'Output directory: '$output_path
echo 'Gwas IDs to analyse:' $*
run_pheno_corr \
	--gwas-ids $(echo $* | tr ' ' ',')\
	--output-path ${output_path}phen_corr_res.txt > ${output_path}01.log
