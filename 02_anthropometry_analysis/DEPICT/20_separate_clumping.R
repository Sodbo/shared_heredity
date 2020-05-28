# Aim of this script is to obtain clumping results separatly for each original trait and for SH. Results are files with rs_ids, using for launch of DEPICT.
# Sript used not unified GWASes and not working properly. It was written to use old version of DEPICT which uses list of rs_id as a start argument. During developing the script it was a decision to not it use further, that's why this script was not finished to working verision.

if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('pROC')) install.packages('pROC'); library('pROC')

path_for_output_results="../../../data/01_anthropometry_results/five_traits/separate_clumping_and_enrichment/"
dir.create(path_for_output_results, showWarnings = FALSE)
#GC corrected original GWASes
gwas_files<-list.files('../../../data/01_anthropometry_results/five_traits/GWAS', full.names=T, pattern='ID_\\d+.csv')
gwas<-lapply(gwas_files, fread)
names(gwas)<-regmatches(gwas_files, regexpr('(?<=ID_)(\\d+)(?=\\.csv)', gwas_files, perl=T))

#GWAS for shared heredity
gwas_sh<-fread('../../../data/01_anthropometry_results/five_traits/linear_combination/SH_GWAS.txt')

gwas<-c(gwas, list(gwas_sh))
names(gwas)[length(gwas)]='SH'

#eaf filtering for minor allele freqency

gwas<-lapply(gwas, function(x){
					tmp=x;
					ind=which(tmp$eaf>=0.01 & (1-tmp$eaf)<=0.99);
					data.frame(tmp[ind,]);
					})

source("../../00_core_functions/clumping.R")
#(locus_table,p_value="p_ma",pos="bp",snp="rs_id", delta=5e5,chr="chr",thr=5e-8,trait=NULL)



out <- lapply(gwas, function(x) function_for_shlop_29_03_2020(locus_table=x,p_value='p',pos='bp',snp='rs_id', delta=5e5, chr='chr', thr=5e-8, trait=NULL))

lapply(1:length(gwas), function(x) fwrite(list(out[[x]]$rs_id),paste0(path_for_output_results,'rs_id_list_for_',names(gwas)[x],'.txt'), col.names=F,row.names=F)
#ind<-lapply(rs_id, function(x) match(snps,x))
#gwas_reordered<-lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]],])

#ind<-match(snps,gwas_sh$SNP)
#gwas_sh=gwas_sh[ind,]

# Forming the input tables 
#pval_orig <- sapply(gwas_reordered, function(x) x$p_gc)
#pval_sh=gwas_sh$p
#chr=gwas_sh$chr
#pos=gwas_sh$pos

#running the main clump_function

#pval_orig=pval_orig;pval_sh=pval_sh;chr=chr;pos=pos;thr_sh_hits=1e-5;thr_sh_sig=5e-8;delta=5e5;SNPs=snps;
#	orig_traits_names=input_file_name;path_output=path_for_output_results



	


