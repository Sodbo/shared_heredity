#this script to reupload the original gwases into test db to
#apply ldscore for calculation of genetic correlations
current_path='../../data/02_Lipids/old_results/GWAS/'
descr_folder=$current_path'01_descriptors'
uni_folder=$current_path'02_unification_out'

for tr in 1287001 1287003 1287004
do
echo Unification of gwas_id=$tr is starting...
mkdir -p ${uni_folder}/ID_$tr

run_uni_qc_rep \
--gwas-path=${current_path}ID_${tr}.csv \
--mapping-path=${current_path}mapping.json \
--descriptors-path=${descr_folder}/descriptor_${tr}.json \
--output-dir=${uni_folder}/ID_${tr}/ \
--output-file=ID_$tr \
--qc-report
#--filter-path filter.txt 

run_upload \
--gwas-path=${uni_folder}/ID_${tr}/ID_${tr}_done.csv

done
