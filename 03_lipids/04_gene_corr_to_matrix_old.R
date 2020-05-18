# Aim of this script to build genetic correlation matrix for traits

path_to_result_directory<-'~/polyomica/projects/shared_heredity/data/02_Lipids/gene_corr/'
path_to_result_directory1<-'~/polyomica/projects/shared_heredity/data/02_Lipids/'

corr_files<-list.files(paste0(path_to_result_directory), full.names=T)
corr_tables<-lapply(corr_files,read.csv)
library(data.table)
united_corr_table<-rbindlist(corr_tables)
united_corr_table$cov<-united_corr_table$rg*sqrt(united_corr_table$h2_obs_1)*sqrt(united_corr_table$h2_obs_2)
united_corr_table$cov_se<-united_corr_table$cov/united_corr_table$z
united_corr_table$cov_se[united_corr_table$gwas_id_1==united_corr_table$gwas_id_2]<-united_corr_table$h2_obs_se_1[united_corr_table$gwas_id_1==united_corr_table$gwas_id_2]
setDT(united_corr_table)
corr_matrix<-dcast(united_corr_table, gwas_id_1 ~ gwas_id_2, value.var = "rg")
colnames(corr_matrix)<-c('',colnames(corr_matrix)[c(-1)])
write.table(corr_matrix, paste0(path_to_result_directory1,'gene_corr_matrix.txt'), row.names=F, quote=F)
cov_matrix<-dcast(united_corr_table, gwas_id_1 ~ gwas_id_2, value.var = "cov")
colnames(cov_matrix)<-c('',colnames(cov_matrix)[c(-1)])
write.table(cov_matrix, paste0(path_to_result_directory1,'gene_cov_matrix.txt'), row.names=F, quote=F)
cov_se_matrix<-dcast(united_corr_table, gwas_id_1 ~ gwas_id_2, value.var = "cov_se")
colnames(cov_se_matrix)<-c('',colnames(cov_se_matrix)[c(-1)])
write.table(cov_se_matrix, paste0(path_to_result_directory1,'gene_cov_se_matrix.txt'), row.names=F, quote=F)
