# Script for unification for GWASes of med Q 
# PGPs done on SOCCS_controls samples

descr_folder='.'
uni_folder=''


for gp in {76..117}
do
echo $gp

mkdir "$uni_folder"/pgp"$gp"

run_uni_qc_rep \
--gwas-path=/mnt/home/upload3/SOCCS.CONTROLS/SOCCS.CONTROLS.PGP"$gp"_A.20190705.txt.gz \
--mapping-path=mapping.json \
--descriptors-path="$descr_folder"/pgp"$gp".json \
--output-dir="$uni_folder"/pgp"$gp"/ \
--output-file=pgp"$gp" \
--filter-path filter.txt \
--qc-report

done
