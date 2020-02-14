descr_folder='/home/ubuntu/polyomica/projects/plasma_v2/03_cohort_gwas_upload/03_SOCCS_controls/01_descriptors'
uni_folder='/home/ubuntu/polyomica/projects/plasma_v2/03_cohort_gwas_upload/03_SOCCS_controls/02_unification_out'

for gp in {1..1}
do
echo $gp

mkdir "$uni_folder"/pgp"$gp"

run_uni_qc_rep \
--gwas-path=~/polyomica/upload3/SOCCS.CONTROLS/SOCCS.CONTROLS.PGP"$gp"_A.20190705.txt.gz \
--mapping-path=mapping.json \
--descriptors-path="$descr_folder"/pgp"$gp".json \
--output-dir="$uni_folder"/pgp"$gp"/ \
--output-file=pgp"$gp" \
--filter-path filter.txt \
--qc-report

done

run_upload \
--gwas-path="$uni_folder"/pgp"$gp"/pgp"$gp"_done.csv


