# genetic covariacne of linear combination with a coeficient with i trait
cov_gi_alpha=function(a,i,covm=covm){
  return(sum(covm[i,]*a))
} 
