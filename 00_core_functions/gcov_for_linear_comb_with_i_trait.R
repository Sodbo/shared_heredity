# genetic covariacne of linear combination with a coeficient with i trait
cov_gi_alfa=function(a,i,covm=covm){
  cov_gi_sum_giai=sum(covm[i,]*a)
  var_gi=covm[i,i]
  var_sum_giai=sum(covm*(a%o%a)) #gen(a)
  
  cor_g_a=cov_gi_sum_giai
  return(cor_g_a)
} 
