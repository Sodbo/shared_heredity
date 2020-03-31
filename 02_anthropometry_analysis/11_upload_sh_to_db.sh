#this script upload SH gwas into test db in order to apply gc_correction

descr_folder='../../data/01_anthropometry_results/linear_combination/01_descriptors'
uni_folder='../../data/01_anthropometry_results/linear_combination/02_unification_out'


mkdir ${uni_folder}/SH

run_uni_qc_rep \
--gwas-path=../../data/01_anthropometry_results/linear_combination/SH_GWAS.txt \
--mapping-path=../../data/01_anthropometry_results/linear_combination/mapping.json \
--descriptors-path=${descr_folder}/descriptor_SH.json \
--output-dir=${uni_folder}/SH/ \
--output-file=SH \
--qc-report
#--filter-path filter.txt 

run_upload \
--gwas-path=${uni_folder}/SH/SH_done.csv

