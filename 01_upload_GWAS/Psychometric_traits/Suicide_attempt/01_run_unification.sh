#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Suicide_attempt/00_raw_data/SA_in_MDD_BIP_SCZ_2019 \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Suicide_attempt/01_upload/sa_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Suicide_attempt/01_upload/sa_descriptor.json \
	--filter-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Suicide_attempt/01_upload/sa_filter.txt \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Suicide_attempt/02_unification_results/ \
	--qc-report \
	--output-file=sa_output

