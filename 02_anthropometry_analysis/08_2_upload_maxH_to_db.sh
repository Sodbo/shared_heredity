#this script upload maxH gwas into test db in order to apply gc_correction
path='../../data/01_anthropometry_results/five_traits/linear_combination/'
descr_folder=$path'01_descriptors'
uni_folder=$path'02_unification_out'


mkdir -p ${uni_folder}/maxH

run_uni_qc_rep \
--gwas-path=${path}maxH_GWAS.txt \
--mapping-path=${path}mapping.json \
--descriptors-path=${descr_folder}/descriptor_SH.json \
--output-dir=${uni_folder}/maxH/ \
--output-file=maxH \
--qc-report
#--filter-path filter.txt 

run_upload \
--gwas-path=${uni_folder}/maxH/maxH_done.csv

