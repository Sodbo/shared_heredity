#gcovm - matrix of genetics covariances
#phem - phneotypic covs
#se - matrix of se of h2 and gcovm
#lll - shrinkage factor [0;1]
#N_permut - N of permutations

source("05_shared_heredity.R", chdir = F)
se <-as.matrix(read.table('~/polyomica/projects/shared_heredity/data/02_Lipids/gene_cov_se_matrix.txt', check.names=F))
#gcovm <-read.table('~/polyomica/projects/shared_heredity/data/02_Lipids/gene_cov_matrix.txt', check.names=F)
#phem <-read.table('~/polyomica/projects/shared_heredity/data/02_Lipids/pheno_corr_matrix.txt', check.names=F)

function_for_estimation_of_alfa_CI=function(gcovm,phem,se,lll=0,N_permut=1000){
  
  #some fucntions
  noise_sd=function(se){
    n=nrow(se)
    out=array(NA,c(n,n))
    for (i in 1:(n)){
      j=i
      for (j in (i):n){
        out[i,j]=rnorm(1,mean=0,sd=se[i,j])
      }
    }
    out[lower.tri(out)]=0
    out_1=t(out)
    out=out+out_1
    diag(out)=diag(out)/2
    out
  }
  
  
  #core
  
  out=array(NA,c(N_permut,nrow(gcovm)))

  i=1
  for (i in 1:N_permut){
    
    NS=noise_sd(se)
    Z_p=gcovm+NS
    
    #shrinkage  
    lambda=rep((1-lll),nrow(Z_p))
    lambda=sqrt(lambda)%o%sqrt(lambda)
    diag(lambda)=1
    Z_p=Z_p*lambda
    
    #alfa estimation!
    alfas=shared_heredity(CovGenTr=Z_p, CorPhenTr=phem)$alphas[2,]
    
    out[i,]=alfas
  }
    
  j=1
  CIs=rep(0,nrow(gcovm))
  for (j in 1:nrow(gcovm)){
    qq=quantile(x=out[,j],probs = c(0.025,0.975))
    CIs[j]=abs(qq[1]-qq[2])/2
    names(CIs)<-colnames(gcovm)
  }
  return(CIs)
}

  (res<-function_for_estimation_of_alfa_CI(A0,CorPhenTr,se,N_permut = 3))
    wanames(res)<-colnames()
write.table(res,'~/polyomica/projects/shared_heredity/data/02_Lipids/CIs_for_3_traits.txt',quote=F)

