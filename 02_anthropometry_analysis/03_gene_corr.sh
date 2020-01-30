path=$1
shift
for i in $*
do
	run_ldscore --rg --gwas-id-1=${i} --gwas-id-2=$(echo $* |tr ' ' ',') --overwrite
done

output_dir=$path'gene_corr'
mkdir $output_dir

for i in $*
do 
	run_ldscore_report --rg --gwas-id-1=${i} --gwas-id-2=$(echo $* |tr ' ' ',') --output-dir $output_dir --output-file gene_corr_$i.txt
done
