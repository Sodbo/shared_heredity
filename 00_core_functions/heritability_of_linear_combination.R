H2=function(a,covm=covm,phem=phem){
  h=sum(covm*(a%o%a))/sum(phem*(a%o%a))
  return(h)
} 
