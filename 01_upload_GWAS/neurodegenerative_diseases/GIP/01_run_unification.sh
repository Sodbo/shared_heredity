#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/four_traits/GIP/00_raw_data/GIP1_GWAS.txt \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/three_traits/linear_combination/01_upload/sh_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/four_traits/GIP/01_upload/gip1_descriptor.json \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/several_traits/four_traits/GIP/02_unification_results/ \
	--qc-report \
	--output-file=gip1_output

