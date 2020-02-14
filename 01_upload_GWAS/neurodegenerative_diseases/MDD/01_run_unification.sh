#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/MDD/00_raw_data/daner_pgc_mdd_NoUKB_No23andMe.txt \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/MDD/01_upload/mdd_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/MDD/01_upload/mdd_descriptor.json \
	--filter-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/MDD/01_upload/mdd_filter.txt \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/MDD/02_unification_results/ \
	--qc-report \
	--output-file=mdd_output

