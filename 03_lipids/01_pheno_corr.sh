function join_by { local IFS=","; echo "$*"; }
run_pheno_corr \
	--gwas-ids $(join_by $*) \
	--output-path ~/polyomica/projects/shared_heredity/data/02_Lipids/phen_corr_res.txt > 01.log
