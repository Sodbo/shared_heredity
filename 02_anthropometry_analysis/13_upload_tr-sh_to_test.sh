#this script upload the gwases "traitX-sh" into test db to
#apply ldscore for calculation of gene correlations between
#gwases without shared heredity component

descr_folder='../../data/anthropometry_results/four_traits/linear_combination/01_descriptors'
uni_folder='../../data/anthropometry_results/four_traits/linear_combination/02_unification_out'

for tr in 4049 4050 4058 4179
do
tr=tr${tr}-SH
echo $tr

mkdir ${uni_folder}/$tr

run_uni_qc_rep \
--gwas-path=../../data/anthropometry_results/four_traits/linear_combination/${tr}_GWAS.txt \
--mapping-path=../../data/anthropometry_results/four_traits/linear_combination/mapping.json \
--descriptors-path=${descr_folder}/descriptor_${tr}.json \
--output-dir=${uni_folder}/${tr}/ \
--output-file=$tr \
#--filter-path filter.txt 
--qc-report

run_upload \
--gwas-path=${uni_folder}/${tr}/${tr}_done.csv

done
