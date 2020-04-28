#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Alcohol/00_raw_data/pgc_alcdep.afr_unrel_genotyped.aug2018_release.txt/pgc_alcdep.afr_unrel_genotyped.aug2018_release.txt \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Alcohol/01_upload/alcohol_afr_unrel_genotyped_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Alcohol/01_upload/alcohol_descriptor.json \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/Alcohol/02_unification_results/ \
	--qc-report \
	--output-file=alcohol_afr_unrel_genotyped_output

