GWAS_linear_combination_Z_based=function(a,Z,covm,N,eaf){
	snp_row<-1:length(Z[,1])
	
	var_y=rep(1,length(a))
	
	se_a<-t(sapply(snp_row, function(x){
								n=N[x,];
								z=Z[x,];
								se_x=sqrt(1/(z^2+n));
								se_x
								} ))
	beta_a<-t(sapply(snp_row, function(x){
								n=N[x,];
								z=Z[x,];
								se_x=se_a[x,];
								b_x=z*se_x
								b_x
								} ))		
	
	out=GWAS_linear_combination_v2(a=a,beta_a=beta_a,se=se_a,var_y=NULL,covm=covm,N=N)
	
	varg=2*(1-eaf)*eaf
	out[,"b"]=out[,"b"]/sqrt(varg)
	out[,"se"]=out[,"se"]/sqrt(varg)
	
	return(out)
} 

GWAS_linear_combination_v2=function(a,beta_a,se,var_y=NULL,covm,N){
	snp_row<-1:length(beta_a[,1])
	if (is.null(var_y)){
		var_y=rep(1,length(a))
	}
	i=1
	for (i in 1:length(var_y)){
		beta_a[,i]=beta_a[,i]/sqrt(var_y[i])
		se[,i]=se[,i]/sqrt(var_y[i])
	}
	
	covm=(1/(sqrt(diag(covm))%o%sqrt(diag(covm))))*covm
		
	#N_min_index=sapply(snp_row, function(x) which.max(N[x,]))
	#se_N_min<-sapply(snp_row, function (x) se[x,N_min_index[x]])
	#N_min<-sapply(snp_row, function(x) N[x,N_min_index[x]])
	#beta_N_min<-sapply(snp_row, function(x) beta_a[x,N_min_index[x]])
	
	#var_y_N_min<-sapply(snp_row, function(x) var_y[N_min_index[x]])
	
	var_y_a=sum(covm*(a%o%a))
	b=(beta_a%*%a)
	varb=sapply(snp_row, function(x){
							sei=se[x,];
							out=sum(covm*(sei%o%sei)*(a%o%a));
							out;
							})
	
	sen=sqrt(varb)
	
	b=b/sqrt(var_y_a)
	sen=sen/sqrt(var_y_a)
	
	N_geom_min=apply(N,MAR=1,function(x){prod(x)^(1/length(x))})
	N_ariph_min=apply(N,MAR=1,function(x){mean(x)})
	N_min=apply(N,MAR=1,function(x){min(x)})
	N_max=apply(N,MAR=1,function(x){max(x)})
	ind=which(abs(b/sen)<2)
	N_from_se=median(1/(sen[ind]^2))
	
	out=cbind(b=b,se=sen,N=N_from_se,N_geom_mean=N_geom_min,N_mean=N_ariph_min,N_min,N_max)
	colnames(out)=c("b","se","N","N_geom_mean","N_mean","N_min","N_max")
	out=as.data.frame(out)
	return(out)
}

