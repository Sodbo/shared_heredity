# genetic covaiance of linear combination with a coefficient a1 and a coefficient a2
cov_gi_a1_a2=function(a1,a2,covm=covm){
  cov_gi_sum_g_a1_a2=sum(covm*(a1%o%a2))
  var_g_a1=sum(covm*(a1%o%a1))
  var_g_a2=sum(covm*(a2%o%a2))
  
  cor_g_a=cov_gi_sum_g_a1_a2
  return(cor_g_a)
} 
