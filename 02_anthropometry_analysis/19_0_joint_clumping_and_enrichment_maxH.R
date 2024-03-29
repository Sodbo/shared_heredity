# Aim of this script is to obtain all figures and tables for MaxH enrichment

if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('pROC')) install.packages('pROC'); library('pROC')
path_to_maxH='../../data/01_anthropometry_results/five_traits/linear_combination/maxH_GWAS_gc_corrected.txt'
path_for_output_results="../../data/01_anthropometry_results/five_traits/joint_clumping_and_enrichment_new/maxH/"
dir.create(path_for_output_results, showWarnings = FALSE)
#GC corrected original GWASes
gwas_files<-list.files('../../data/01_anthropometry_results/five_traits/GWAS', full.names=T, pattern='ID_\\d+_gc_corrected.csv')

#The last gwas is for shared heredity and will not be used
gwas_files<-gwas_files[-6]

#GWASes for original traits
gwas<-lapply(gwas_files, fread)
names(gwas)<-regmatches(gwas_files, regexpr('(?<=ID_)(\\d+)(?=_gc_corrected\\.csv)', gwas_files, perl=T))

#GWAS for maxH
gwas_maxH=fread(path_to_maxH)



#eaf filtering for 
ind=which(gwas_maxH$eaf>=0.01 & gwas_maxH$eaf<=0.99)
gwas_maxH=gwas_maxH[ind,]

gwas<-lapply(gwas, function(x){
					tmp=x;
					ind=which(tmp$eaf>=0.01 & tmp$eaf<=0.99);
					tmp[ind,];
					} )


#reordering of all original GWASs and maxH
rs_id<-lapply(gwas, function(x) x$rs_id)
snps<-Reduce(intersect,rs_id)
snps<-intersect(snps,gwas_maxH$SNP)

ind<-lapply(rs_id, function(x) match(snps,x))
gwas_reordered<-lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]],])

ind<-match(snps,gwas_maxH$SNP)
gwas_maxH=gwas_maxH[ind,]

# Forming the input tables 
pval_orig <- sapply(gwas_reordered, function(x) x$p_gc)
pval_sh=gwas_maxH$p
chr=gwas_maxH$chr
pos=gwas_maxH$pos

#running the main clump_function

#pval_orig=pval_orig;pval_sh=pval_sh;chr=chr;pos=pos;thr_sh_hits=1e-5;thr_sh_sig=5e-8;delta=5e5;SNPs=snps;
#	orig_traits_names=input_file_name;path_output=path_for_output_results


source("../00_core_functions/joint_function_for_enrichment_and_auc.R")
	
out=joint_function(pval_orig=pval_orig,pval_sh=pval_sh,chr=chr,pos=pos,thr_sh_hits=1e-5,thr_sh_sig=5e-8,
	orig_traits_names=names(gwas),path_output=path_for_output_results,delta=5e5,SNPs=snps)
