path=$1;
export PROD=T
run_ldscore --h2 --gwas-id=4049,4050,4058,4179

\copy (select * from gwas.ldsc where gwas_id in (4049,4050,4058,4179);) To ${path}intrecept_data.csv' With CSV HEADER DELIMITER ',';
