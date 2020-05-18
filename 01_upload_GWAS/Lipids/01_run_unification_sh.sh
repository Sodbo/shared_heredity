#!/bin/bash

run_uni_qc_rep \
	--gwas-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/00_raw_data/SH_Lipids_GWAS.txt \
	--mapping-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/01_upload/sh_mapping.json \
	--descriptors-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/01_upload/sh_descriptor.json \
	--output-dir=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/02_unification_results/ \
	--qc-report \
	--output-file=sh_output

