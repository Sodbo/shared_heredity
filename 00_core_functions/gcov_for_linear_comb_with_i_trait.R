# genetic covariance of linear combination with a coefficient "a" with trait "i"
cov_gi_alpha=function(a,i,covm=covm){
  return(sum(covm[i,]*a))
} 
