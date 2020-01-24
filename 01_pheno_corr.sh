function join_by { local IFS=","; echo "$*"; }
run_pheno_corr \
	--gwas-ids $(join_by $*) \
	--output-path ~/polyomica/projects/plasma_v2/tiys_src/shared_heredity/phen_corr_res.txt > 01.log
