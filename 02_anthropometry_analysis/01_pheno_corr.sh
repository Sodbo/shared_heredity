output_path=$1
shift
echo 'Output directory: '$output_path
echo 'Gwas IDs to analyse:' $*
python3 ~/code_folder/gwas_v2/gwas_v2/analysis/database_interfaces/run_phen_corr.py \
	--gwas-ids $(echo $* | tr ' ' ',')\
	--output-path ${output_path}phen_corr_res.txt > ${output_path}01.log
