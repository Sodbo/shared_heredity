#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/00_raw_data/happiness.csv \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/01_upload/happiness_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/01_upload/happiness_descriptor.json \
	--filter-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/01_upload/happiness_filter.txt \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/happiness/02_unification_results/ \
	--qc-report \
	--output-file=happiness_output

