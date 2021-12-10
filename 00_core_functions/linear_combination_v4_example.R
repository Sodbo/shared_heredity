library(Rcpp)
sourceCpp("./linear_combination_v4.cpp")

linear_combination_v4 <- function(a, Z, covm, N, eaf) {    

  if(class(a) == "numeric") { # a = vector
    result1 <- GWAS_linear_combination_Z_basedC(a, Z, covm, N, eaf)
    return(result1)

  } else { # a = matrix
  
    results <- GWAS_linear_combination_Z_basedC_GIPs(a, Z, covm, N, eaf)
    result_list <- list()
    
    for(i in 1:ncol(a)) {
      result_list[[paste0("GIP",i)]] <- results[,(1+(i-1)*3) : (i*3)]
      names(result_list[[paste0("GIP",i)]]) <- c("b","se","N")
    }
    
    return(result_list)  

  }

}

## START

load("st03_01_loadings.Rdata") # 'res' variable

pcov_mat <- res$GIPs$cov_y
gips <- res$GIPs$GIP_coeff

gip_stats <- linear_combination_v4(a=gips, Z, covm, N, eaf)
