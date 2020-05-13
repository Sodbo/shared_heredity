#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/05_minus_sh/happiness-sh_gwas.txt \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/06_upload_minus_sh/happiness-sh_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/06_upload_minus_sh/happiness-sh_descriptor.json \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/07_minus_sh_unified/ \
	--qc-report \
	--output-file=happiness-sh_output

