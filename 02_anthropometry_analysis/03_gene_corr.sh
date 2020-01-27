function join_by { local IFS=","; echo "$*"; }
for i in $*
do
	run_ldscore --rg --gwas-id-1=${i} --gwas-id-2=$(join_by $*) --overwrite
done

output_dir='gene_corr'
mkdir $output_dir

for i in $*
do 
	run_ldscore_report --rg --gwas-id-1=${i} --gwas-id-2=$(join_by $*) --output-dir $output_dir --output-file gene_corr_$i.txt
done
