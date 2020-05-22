#!/bin/bash

run_uni_qc_rep \
	--gwas-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/GIP/00_raw_data/GIP1_GWAS.txt \
	--mapping-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/01_upload/sh_mapping.json \
	--descriptors-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/GIP/01_upload/gip1_descriptor.json \
	--output-dir=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/GIP/02_unification_results/ \
	--qc-report \
	--output-file=gip1_output

