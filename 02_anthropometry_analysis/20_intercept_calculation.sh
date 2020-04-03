{
	output_path=../../data/01_anthropometry_results/Traits_minus_SH/
source ../00_core_functions/p_correction_for_gc.sh $output_path 181 182 183 184 185
} 2>&1 | tee $output_path/log_20.txt
