library(data.table)
library(dplyr)

path_to_result_directory<-'../data/anthropometry_results/four_traits/GWAS/scaled_filtered/'
gwas_files<-list.files(path_to_result_directory, full.names=T, pattern="ID_\\d+.csv")
gwas<-lapply(gwas_files, fread)

aa<-read.table('../data/anthropometry_results/four_traits/alphas.txt', row.names=1)
covm<-read.table('../data/anthropometry_results/four_traits/pheno_corr_matrix.txt', row.names=1,  check.names=F)

rs_id<-lapply(gwas, function(x) x$rs_id)
snps<-Reduce(intersect,rs_id)
ind<-lapply(rs_id, function(x) match(snps,x))

#Is it necessary to remove NA for the case of not whole intersection?
#ind<-lapply(ind,function(x) )

n<-as.list(c(1:length(ind)))

gwas_reordered<-lapply(n, function(x) gwas[[x]][ind[[x]],])

betas<-lapply(gwas_reordered, function(x) x$beta)
betas <- do.call(cbind, betas)

se1 <- gwas_reordered[[1]]$se

GWAS_linear_combination=function(a,betaa,se1,vary1=1,covm,N){
    vary2=sum(covm*(a%o%a))
    b=(betaa%*%a)
    
    vary1_varg_n=se1^2+(betaa[,1]^2)/N
    
    se2=vary1_varg_n*(vary2/vary1)
    
    se2=se2-b^2/N
    
    se=sqrt(se2)
    
    b=b/sqrt(vary2)
    se=se/sqrt(vary2)
    
    out=cbind(b=b,se=se)
    colnames(out)=c("b","se")
    out=as.data.frame(out)
    return(out)
}

x<-lapply(n,

for (NPC in 1:4){
  x=GWAS_linear_combination(a=as.vector(aa[,NPC]),betaa=as.matrix(betas),se1=as.vector(se1),vary1=1,covm=as.matrix(covm),N=336107)
  x=mutate(x,Z=b/se,p=pchisq(Z^2,1,low=F))
  x=mutate(x,SNP=back$rs_id)
  x=mutate(x,A1=back$ea,A2=back$ra,N=336107,chr=back$chr,pos=back$bp,
             eaf=back$eaf)
  
  head(x,n=2)
  
  fnme=paste0("/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/GPC/20181011_DISC_GPC",NPC,".txt")
    
  data.table::fwrite(
    x, 
    row.names=F,
    file = fnme,
    sep = '\t')
  }
  
