path<-commandArgs(trailingOnly=T)[1]
x<-read.csv(paste0(path,'phen_corr_res.txt'))
x_mirrored=x[,c(1,3,2,4,5)]
colnames(x_mirrored)[c(2,3)]=colnames(x_mirrored)[c(3,2)]
x=rbind(x,x_mirrored)
library(data.table)
setDT(x)
corr_matrix=dcast(x, gwas_id_1 ~ gwas_id_2, value.var = "corr")
corr_matrix[is.na(corr_matrix)]=1
colnames(corr_matrix)=c('',colnames(corr_matrix)[c(-1)])
write.table(corr_matrix,paste0(path,'pheno_corr_matrix.txt'), row.names=F, quote=F)
