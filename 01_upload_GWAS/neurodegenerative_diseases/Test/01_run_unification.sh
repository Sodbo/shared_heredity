#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/autism_spectrum_disorder/00_raw_data/daner_AUT_meta14_CEU_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.txt \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/autism_spectrum_disorder/01_upload/asd_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/autism_spectrum_disorder/01_upload/asd_descriptor.json \
	--filter-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/autism_spectrum_disorder/01_upload/asd_filter.txt \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/src/01_upload_GWAS/neurodegenerative_diseases/Test/ \
	--qc-report \
	--output-file=asd_output

