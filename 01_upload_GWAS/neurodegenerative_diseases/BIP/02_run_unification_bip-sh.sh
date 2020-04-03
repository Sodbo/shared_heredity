#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/BIP/05_minus_sh/bip-sh_gwas.txt \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/BIP/06_upload_minus_sh/bip-sh_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/BIP/06_upload_minus_sh/bip-sh_descriptor.json \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/BIP/07_minus_sh_unified/ \
	--qc-report \
	--output-file=bip-sh_output

