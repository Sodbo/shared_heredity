GWAS_linear_combination_v3=function(a,Z,covm,N){
	snp_row<-1:length(Z[,1])
	
	var_y=rep(1,length(a))
	
	se_a<-t(sapply(snp_row, function(x){
								n=N[x,];
								z=Z[x,];
								se_x=sqrt(1/(z^2+n^2));
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
	return(out)
} 
