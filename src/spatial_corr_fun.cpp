#include "RcppArmadillo.h"
#include "Spillover.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List spatial_corr_fun(double phi,
                            arma::mat spatial_dists){

int p_z = spatial_dists.n_cols;
double log_deter = 0.00; 
double sign = 0.00;     
arma::mat spatial_corr(p_z, p_z); spatial_corr.fill(0);

for(int j = 0; j < p_z; ++j){
   for(int k = 0; k < p_z; ++k){
      if(spatial_dists(j,k) <= (1.00/phi)){
        spatial_corr(j,k) = 1.00 - 
                            1.50*phi*spatial_dists(j,k) + 
                            0.50*(phi*spatial_dists(j,k))*(phi*spatial_dists(j,k))*(phi*spatial_dists(j,k));
        }
      }
   }

arma::mat spatial_corr_inv = inv_sympd(spatial_corr);
log_det(log_deter, sign, spatial_corr);

return Rcpp::List::create(Rcpp::Named("spatial_corr_inv") = spatial_corr_inv,
                          Rcpp::Named("log_deter") = log_deter);

}

