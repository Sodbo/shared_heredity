#this script upload the gwases "traitX-sh" into test db to
#apply ldscore for calculation of genetic correlations between
#gwases without shared heredity component
path='../../data/01_anthropometry_results/five_traits/linear_combination/'
descr_folder=$path'01_descriptors'
uni_folder=$path'02_unification_out'

for tr in 191 192 193 194 199
do
tr=Tr${tr}-SH
echo $tr

mkdir ${uni_folder}/$tr

run_uni_qc_rep \
--gwas-path=$path${tr}_GWAS.txt \
--mapping-path=${path}mapping.json \
--descriptors-path=${descr_folder}/descriptor_${tr}.json \
--output-dir=${uni_folder}/${tr}/ \
--output-file=$tr \
--qc-report
#--filter-path filter.txt 

run_upload \
--gwas-path=${uni_folder}/${tr}/${tr}_done.csv

done
