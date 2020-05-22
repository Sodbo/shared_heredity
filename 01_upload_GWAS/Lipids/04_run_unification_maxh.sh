#!/bin/bash

run_uni_qc_rep \
	--gwas-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/MaxH/00_raw_data/MaxH_GWAS.txt \
	--mapping-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/SH/01_upload/sh_mapping.json \
	--descriptors-path=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/MaxH/01_upload/maxh_descriptor.json \
	--output-dir=/mnt/polyomica/projects/shared_heredity/data/02_Lipids/three_traits/MaxH/02_unification_results/ \
	--qc-report \
	--output-file=maxh_output

