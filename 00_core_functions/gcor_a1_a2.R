# genetic correlation of linear combination with a coeeficient a1 and a coeeficient a2
cor_gi_a1_a2=function(a1,a2,covm=covm){
  cov_gi_sum_g_a1_a2=sum(covm*(a1%o%a2))
  var_g_a1=sum(covm*(a1%o%a1))
  var_g_a2=sum(covm*(a2%o%a2))
  
  cor_g_a=cov_gi_sum_g_a1_a2/sqrt(var_g_a1*var_g_a2)
  return(cor_g_a)
}

 
