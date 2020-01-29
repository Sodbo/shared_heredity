output_path=$1
shift
function join_by_comma { local IFS=","; echo $*; }
run_pheno_corr \
	--gwas-ids $(join_by_comma $*) \
	--output-path ${output_path}/phen_corr_res.txt > 01.log
