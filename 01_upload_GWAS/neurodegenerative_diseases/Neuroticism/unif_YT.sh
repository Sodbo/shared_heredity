
WORKDIR=/home/ubuntu/neurotism_genable

run_uni_qc_rep \
--gwas-path=/home/ubuntu/polyomica/projects/neurotism_genabase/sumstats_neuro_sum_ctg_format.txt \
--mapping-path=/home/ubuntu/polyomica/projects/neurotism_genabase/neu_mapping_YT.json \
--descriptors-path=/home/ubuntu/polyomica/projects/neurotism_genabase/neu_descriptor_YT_n.json \
--output-dir=/home/ubuntu/polyomica/projects/shared_heredity/data/03_neurodegenerative_diseases/Neuroticism/01_unification_results/ \
--output-file=neurot \
--filter-path=/home/ubuntu/polyomica/projects/neurotism_genabase/filter.txt \
--qc-report 

#run_unifier --gwas-path=sumstats_neuro_sum_ctg_format.txt --mapping-path=neu_mapping.json --descriptors-path=neu_descriptor.json --output-dir=/home/ubuntu/neurotism_genable

