GWAS_linear_combination=function(a,beta_a,se,var_y=rep(1,length(a)),covm,N){
	snp_row<-1:length(beta_a[,1])
	
	i=1
	for (i in 1:length(var_y)){
		beta_a[,i]=beta_a[,i]/sqrt(var_y[i])
		se[,i]=se[,i]/sqrt(var_y[i])
	}
	
	covm=(1/(sqrt(diag(covm))%o%sqrt(diag(covm))))*covm
	
	N_min_index=sapply(snp_row, function(x) which.min(N[x,]))
	se_N_min<-sapply(snp_row, function (x) se[x,N_min_index[x]])
	N_min<-sapply(snp_row, function(x) N[x,N_min_index[x]])
	beta_N_min<-sapply(snp_row, function(x) beta_a[x,N_min_index[x]])
	
	#var_y_N_min<-sapply(snp_row, function(x) var_y[N_min_index[x]])
	
	var_y_a=sum(covm*(a%o%a))
	b=(beta_a%*%a)
	
	var_y_N_min_to_var_g_N_min_ratio=se_N_min^2+(beta_N_min^2)/N_min
	
	#se2=var_y_N_min_to_var_g_N_min_ratio*(var_y_a/var_y_N_min)
	se2=var_y_N_min_to_var_g_N_min_ratio*(var_y_a/1)
	
	se2=se2-b^2/N_min
	
	se=sqrt(se2)
	
	b=b/sqrt(var_y_a)
	se=se/sqrt(var_y_a)
	
	out=cbind(b=b,se=se,N=N_min)
	colnames(out)=c("b","se","N")
	out=as.data.frame(out)
	return(out)
}
