#!/bin/bash

python3 /home/ubuntu/code_folder/gwas_v2/gwas_v2/unifier/run_uni_qc_rep.py \
	--gwas-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/BIP/00_raw_data/daner_PGC_BIP32b_mds7a_0416a.txt \
	--mapping-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/BIP/01_upload/bip_mapping.json \
	--descriptors-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/BIP/01_upload/bip_descriptor.json \
	--filter-path=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/BIP/01_upload/bip_filter.txt \
	--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/neurodegenerative_diseases/BIP/02_unification_results/ \
	--qc-report \
	--output-file=bip_output

