#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/SCZ/00_raw_data/ckqny.scz2snpres \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/SCZ/01_upload/scz_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/SCZ/01_upload/scz_descriptor.json \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/SCZ/02_unification_results/ \
	--qc-report \
	--output-file=scz_output

