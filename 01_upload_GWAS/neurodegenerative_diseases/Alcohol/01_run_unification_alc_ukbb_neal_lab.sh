#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/alcohol/00_raw_data/Alcohol_intake_frequency.csv \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/alcohol/01_upload/alc_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/alcohol/01_upload/alc_descriptor.json \
	--filter-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/alcohol/01_upload/alc_filter.txt \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/alcohol/02_unification_results/ \
	--qc-report \
	--output-file=alc_output

