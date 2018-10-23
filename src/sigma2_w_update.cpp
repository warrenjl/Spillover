#include "RcppArmadillo.h"
#include "Spillover.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_w_update(arma::vec w,
                       arma::mat corr_inv,
                       double alpha_sigma2_w,
                       double beta_sigma2_w){

int p_w = w.size();
  
double alpha_sigma2_w_update = 0.50*p_w + 
                               alpha_sigma2_w;

double beta_sigma2_w_update = 0.50*dot(w, ((corr_inv)*w)) + 
                              beta_sigma2_w;

double sigma2_w = 1/R::rgamma(alpha_sigma2_w_update,
                              (1/beta_sigma2_w_update));

return(sigma2_w);

}





