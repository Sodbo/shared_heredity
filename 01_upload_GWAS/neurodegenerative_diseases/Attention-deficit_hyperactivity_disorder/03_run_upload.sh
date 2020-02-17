# Script to upload GWASes of med Q PGPs 

descr_folder='/home/ubuntu/polyomica/projects/plasma_v2/03_cohort_gwas_upload/03_SOCCS_controls/01_descriptors'
uni_folder='/home/ubuntu/polyomica/projects/plasma_v2/03_cohort_gwas_upload/03_SOCCS_controls/02_unification_out'

for gp in {76..117}
do
echo $gp

run_upload \
--gwas-path="$uni_folder"/pgp"$gp"/pgp"$gp"_done.csv

done

