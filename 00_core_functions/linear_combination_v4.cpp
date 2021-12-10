#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::DataFrame GWAS_linear_combination_Z_basedC(NumericVector a, NumericMatrix Z, NumericMatrix covm, NumericMatrix N, NumericVector eaf)
{
  
  int p = 0, progress = 0;
  int nrow = Z.nrow(), ncol = Z.ncol();
  NumericMatrix beta(nrow, ncol), se(nrow, ncol);
  NumericVector shared_beta(nrow), shared_se(nrow);
  NumericVector N_geom_min(nrow,1.0), N_ariph_min(nrow,0.0), N_min(nrow,99999999.9), N_max(nrow, 0.0);
  std::vector<double> N_est; // vector to contain estimated N values from null SNPs

  Rcpp::Rcout << "Combining traits... " << std::endl;

  for (int i = 0; i < nrow; i++) { // i = SNP
    
    for (int j = 0; j < ncol; j++) { // j = TRAIT

      se(i,j) = sqrt(1 / ( pow(Z(i,j), 2) + N(i,j) ) ); // standardised stderrs
      beta(i,j) = Z(i,j) * se(i,j); // standardised betas
      
      shared_beta(i) += beta(i,j) * a(j); // combined beta
      shared_se(i) += pow(se(i,j),2) * pow(a(j),2); // unadjusted combined variance

      /// SAMPLE SIZES

      N_geom_min(i) *= N(i,j);
      N_ariph_min(i) += N(i,j)/ncol; 
      if(N(i,j) < N_min(i)) N_min(i) = N(i,j);
      if(N(i,j) > N_max(i)) N_max(i) = N(i,j);

    }

    N_geom_min(i) = pow(N_geom_min(i), 1.0/ncol);


      for (int j1 = 0; j1 < (ncol-1); j1++) {
        for (int j2 = j1+1; j2 < ncol; j2++)
        {
          
          shared_se(i) += 2.0 * a(j1) * a(j2) * se(i,j1) * se(i,j2) * covm(j1,j2); // adjust combined variance
        
        }
      }
    
    shared_se(i) = sqrt(shared_se(i)); // convert to stderr
    
    // Estimate N if null SNP
    double z = shared_beta(i) / shared_se(i); 
    if(abs(z) < 1.959964) N_est.push_back(1/pow(shared_se(i),2)); 


    // adjust for allele frequency
    double varg = sqrt(2.0 * (1.0-eaf(i)) * eaf(i)); 
    shared_beta(i) /= varg;
    shared_se(i) /= varg;

    progress = 100 * (float)i / (float)(nrow);
    if (progress > p)
    {
      p = progress;
      Rcpp::Rcout << "\r" << progress << "%";
    }
    
  }

  Rcpp::Rcout << "\r100%" << std::endl;

  // Get median estimate of N
  sort(N_est.begin(), N_est.end());
  int s = N_est.size();
  double N_median;
  
  if (s % 2 == 0) {
      N_median = (N_est[s/2 - 1] + N_est[s/2]) / 2;
  } else {
      N_median = N_est[s/2];
  }

  return DataFrame::create(
    Named("b") = shared_beta,
    Named("se") = shared_se,
    Named("N") = N_median,
    Named("N_geom_mean") = N_geom_min,
    Named("N_mean") = N_ariph_min,
    Named("N_min") = N_min,
    Named("N_max") = N_max);
}
    








// [[Rcpp::export]]
Rcpp::DataFrame GWAS_linear_combination_Z_basedC_GIPs(NumericMatrix a, NumericMatrix Z, NumericMatrix covm, NumericMatrix N, NumericVector eaf)
{
  int nrow = Z.nrow(), ncol = Z.ncol();
  List results(ncol);

//// STEP 1
  Rcpp::Rcout << "Standardising statistics... " << std::endl;

  NumericMatrix beta(nrow, ncol), se(nrow, ncol);
  int p = 0, progress = 0;

  for (int i = 0; i < nrow; i++) { // i = SNP
      
    for (int j = 0; j < ncol; j++) { // j = TRAIT

    se(i,j) = sqrt(1 / ( pow(Z(i,j), 2) + N(i,j) ) ); // standardised stderrs
    beta(i,j) = Z(i,j) * se(i,j); // standardised betas
    
    }
  
    // report progress
    progress = 100 * (float)i / (float)(nrow);
    if (progress > p)
      {
        p = progress;
        Rcpp::Rcout << "\r" << progress << "%";
      }
  
  }

  Rcpp::Rcout << "\r" << "100%" << std::endl << std::endl;
  

//// STEP 2
  Rcpp::Rcout << "Combining traits... " << std::endl;

  for (int h = 0; h < ncol; h++) { // h = GIP

  int p = 0, progress = 0;
  NumericVector shared_beta(nrow), shared_se(nrow);

  std::vector<double> N_est; // vector to contain estimated N values from null SNPs


    for (int i = 0; i < nrow; i++) { // i = SNP
      
      for (int j = 0; j < ncol; j++) { // j = TRAIT

        shared_beta(i) += beta(i,j) * a(j,h); // combined beta
        shared_se(i) += pow(se(i,j),2) * pow(a(j,h),2); // unadjusted combined variance

      }

        for (int j1 = 0; j1 < (ncol-1); j1++) {
          for (int j2 = j1+1; j2 < ncol; j2++)
          {
            
            shared_se(i) += 2.0 * a(j1,h) * a(j2,h) * se(i,j1) * se(i,j2) * covm(j1,j2); // adjust combined variance
          
          }
        }
      
      shared_se(i) = sqrt(shared_se(i)); // convert to stderr
      
      // Estimate N if null SNP
      double z = shared_beta(i) / shared_se(i); 
      if(abs(z) < 1.959964) N_est.push_back(1/pow(shared_se(i),2)); 


      // adjust for allele frequency
      double varg = sqrt(2.0 * (1.0-eaf(i)) * eaf(i)); 
      shared_beta(i) /= varg;
      shared_se(i) /= varg;
      
      progress = 100 * (float)i / (float)(nrow);
      if (progress > p)
      {
        p = progress;
        Rcpp::Rcout << "\rGIP" << h+1 << " " << progress << "%";
      }
      
    }
    
    Rcpp::Rcout << "\rGIP" << h+1 << " 100%" << std::endl;


    // Get median estimate of N
    sort(N_est.begin(), N_est.end());
    int s = N_est.size();
    double N_median;
    
    if (s % 2 == 0) {
        N_median = (N_est[s/2 - 1] + N_est[s/2]) / 2;
    } else {
        N_median = N_est[s/2];
    }

    results(h) = DataFrame::create(
      Named("b") = shared_beta,
      Named("se") = shared_se,
      Named("N") = N_median);
  }

return results;
}



