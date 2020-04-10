#this script to reupload the original gwases into test db to
#apply ldscore for calculation of genetic correlations

descr_folder='../../data/01_anthropometry_results/GWAS/scaled_filtered/reupload_gwases/01_descriptors'
uni_folder='../../data/01_anthropometry_results/GWAS/scaled_filtered/reupload_gwases/02_unification_out'

for tr in 4049 4050 4058 4179
do
echo Unification of gwas_id=$tr is starting...
mkdir -p ${uni_folder}/ID_$tr

run_uni_qc_rep \
--gwas-path=../../data/01_anthropometry_results/GWAS/scaled_filtered/ID_${tr}.csv \
--mapping-path=../../data/01_anthropometry_results/GWAS/scaled_filtered/reupload_gwases/mapping.json \
--descriptors-path=${descr_folder}/descriptor_${tr}.json \
--output-dir=${uni_folder}/ID_${tr}/ \
--output-file=ID_$tr \
--qc-report
#--filter-path filter.txt 

run_upload \
--gwas-path=${uni_folder}/ID_${tr}/ID_${tr}_done.csv

done
