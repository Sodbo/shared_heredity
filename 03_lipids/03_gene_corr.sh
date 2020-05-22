# This script is to calculate genetic correlatioins between traits
# using tools incorporated into GWAS-Map database

# This script uses command arguments set in 00a_start.sh or 16a_start.sh script
path=$1
shift
output_dir=$path'gene_corr'
#mkdir $output_dir
#for i in $*
#do
	#run_ldscore --rg --gwas-id-1=$i --gwas-id-2=$(echo $* |tr ' ' ',') --overwrite
	#run_ldscore --h2 --gwas-id=$i --overwrite
	#run_ldscore_report --rg --gwas-id-1=$i --gwas-id-2=$(echo $* |tr ' ' ',') --output-dir $output_dir --output-file gene_corr_$i.txt
#done

run_ldscore_report --h2 --gwas-id=$(echo $* |tr ' ' ',') --output-dir $output_dir --output-file h2.txt
