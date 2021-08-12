# Aim of this script is to perform joint clumping on original traits, SGCT and UGCTs

if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('pROC')) install.packages('pROC'); library('pROC')

setwd('/mnt/polyomica/projects/shared_heredity/')
path_for_output_results="./data/01_anthropometry_results/five_traits/joint_clumping_and_enrichment_full_set/"
dir.create(path_for_output_results, showWarnings = FALSE)
# GC corrected GWASes on original traits and SGCT
gwas_files<-list.files('./data/01_anthropometry_results/five_traits/GWAS/', full.names=T, pattern='ID_\\d+_gc_corrected.csv')
# add non-corrected GWASes on UGCTs
ugct_path<-'/mnt/polyomica/projects/shared_heredity/data/01_anthropometry_results/five_traits/linear_combination/02_unification_out/'
ugct_files<-sapply(c(191:194,199), function(x) paste0(ugct_path, 'Tr',x,'-SH/Tr',x,'-SH_done.csv'))
gwas_files<-c(gwas_files,ugct_files)
gwas<-lapply(gwas_files, fread)
names(gwas)<-c('BMI', 'Weight', 'Hip', 'Waist', 'Fat', 'SGCT', 'BMI UGCT', 'Weight UGCT', 'Hip UGCT', 'Waist UGCT', 'Fat UGCT')

#GWAS for SGCT
gwas_sh=gwas$SGCT

#original traits and UGCTs
gwas=gwas[-6]


#eaf filtering for 
ind=which(gwas_sh$eaf>=0.01 & gwas_sh$eaf<=0.99)
gwas_sh=gwas_sh[ind,]

gwas<-lapply(gwas, function(x){
					tmp=x;
					ind=which(tmp$eaf>=0.01 & tmp$eaf<=0.99);
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
pval_orig <- sapply(gwas_reordered[1:5], function(x) x$p_gc)
pval_ugct <- sapply(gwas_reordered[6:10], function(x) x$p)
pval_orig <- cbind(pval_orig, pval_ugct)
pval_sh=gwas_sh$p_gc
chr=gwas_sh$chr
pos=gwas_sh$pos

#running the main clump_function

#pval_orig=pval_orig;pval_sh=pval_sh;chr=chr;pos=pos;thr_sh_hits=1e-5;thr_sh_sig=5e-8;delta=5e5;SNPs=snps;
#	orig_traits_names=input_file_name;path_output=path_for_output_results


source("./elgaeva_src/shared_heredity/00_core_functions/joint_function_for_enrichment_and_auc.R")
	
out=joint_function(pval_orig=pval_orig,pval_sh=pval_sh,chr=chr,pos=pos,thr_sh_hits=1e-5,thr_sh_sig=5e-8,
	orig_traits_names=names(gwas),path_output=path_for_output_results,delta=5e5,SNPs=snps)
