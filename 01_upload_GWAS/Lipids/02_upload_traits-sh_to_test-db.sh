# This script to upload traits minus shared heredity to test db to

current_path='/mnt/polyomica/projects/shared_heredity/data/02_Lipids/traits_minus_SH/GWAS'

for tr in 1287001 1287003 1287004
do
echo Unification of gwas_id=$tr is starting...

run_uni_qc_rep \
--gwas-path=${current_path}/ID_${tr}/00_raw_trait-sh/ID_${tr}-sh_gwas.txt \
--mapping-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/01_upload/sh_mapping.json \
--descriptors-path=${current_path}/ID_${tr}/01_upload_trait-sh/descriptor_${tr}-sh.json \
--output-dir=${current_path}/ID_${tr}/02_unification_results_trait-sh/ \
--output-file=ID_${tr}-sh \
--qc-report
#--filter-path filter.txt 

run_upload \
--gwas-path=${current_path}/ID_${tr}/02_unification_results_trait-sh/ID_${tr}-sh_done.csv

done
