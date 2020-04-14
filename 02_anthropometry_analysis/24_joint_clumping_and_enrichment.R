# Aim of this script is to obtain all figures and tables for SH enrichment

if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('pROC')) install.packages('pROC'); library('pROC')

path_for_output_results="../../data/01_anthropometry_results/GWAS/scaled_filtered/joint_clumping_and_enrichment_new/"
dir.create(path_for_output_results, showWarnings = FALSE)
path_to_gwas='../../data/01_anthropometry_results/GWAS/scaled_filtered/'
#GC corrected original GWASes
input_file_name<-c('ID_4049_gc_corrected.csv','ID_4050_gc_corrected.csv','ID_4058_gc_corrected.csv','ID_4179_gc_corrected.csv')
gwas<-lapply(paste0(path_to_gwas, input_file_name), fread)

#GWAS for shared heredity
gwas_sh=fread("../../data/01_anthropometry_results/linear_combination/SH_GWAS.txt",data.table=F)

#eaf filtering for 
ind=which(gwas_sh$eaf>=0.01 & (1-gwas_sh$eaf)<=0.99)
gwas_sh=gwas_sh[ind,]

gwas<-lapply(gwas, function(x){
					tmp=x;
					ind=which(tmp$eaf>=0.01 & (1-tmp$eaf)<=0.99);
					tmp[ind,];
					} )


#reordering of all original GWASs and SH
rs_id<-lapply(gwas, function(x) x$rs_id)
snps<-Reduce(intersect,rs_id)
snps<-intersect(snps,gwas_sh$SNP)

ind<-lapply(rs_id, function(x) match(snps,x))
gwas_reordered<-lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]],])

ind<-match(snps,gwas_sh$SNP)
gwas_sh=gwas_sh[ind,]

# Forming the input tables 
pval_orig <- sapply(gwas_reordered, function(x) x$p_gc)
pval_sh=gwas_sh$p
chr=gwas_sh$chr
pos=gwas_sh$pos

#running the main clump_function

#pval_orig=pval_orig;pval_sh=pval_sh;chr=chr;pos=pos;thr_sh_hits=1e-5;thr_sh_sig=5e-8;delta=5e5;SNPs=snps;
#	orig_traits_names=input_file_name;path_output=path_for_output_results


source("../00_core_functions/joint_function_for_enrichment_and_auc.R")
	
out=joint_function(pval_orig=pval_orig,pval_sh=pval_sh,chr=chr,pos=pos,thr_sh_hits=1e-5,thr_sh_sig=5e-8,
	orig_traits_names=input_file_name,path_output=path_for_output_results,delta=5e5,SNPs=snps)
