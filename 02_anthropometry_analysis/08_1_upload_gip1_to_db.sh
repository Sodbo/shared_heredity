#this script upload GIP1 gwas into test db in order to apply gc_correction
path='../../data/01_anthropometry_results/five_traits/GIP1/'
descr_folder=$path'01_descriptors'
uni_folder=$path'02_unification_out'


mkdir -p ${uni_folder}/GIP1

run_uni_qc_rep \
--gwas-path=${path}GIP1_GWAS.txt \
--mapping-path=${path}mapping.json \
--descriptors-path=${descr_folder}/descriptor_GIP1.json \
--output-dir=${uni_folder}/GIP1/ \
--output-file=GIP1 \
--qc-report
#--filter-path filter.txt 

run_upload \
--gwas-path=${uni_folder}/GIP1/GIP1_done.csv

