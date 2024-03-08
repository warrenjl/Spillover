#include "RcppArmadillo.h"
#include "Spillover.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_aux_update(arma::vec y,
                        arma::mat x,
                        arma::vec tri_als,
                        arma::vec spillover_covar,
                        arma::mat z,
                        arma::vec beta_old,
                        double lambda_old,
                        arma::vec w_old){
  
arma::vec mean_w_aux = x*beta_old +
                       spillover_covar*lambda_old +
                       z*w_old;

arma::vec input = tri_als;
arma::vec w_aux = rcpp_pgdraw(input,
                              mean_w_aux);

arma::vec gamma = (y - 0.50*tri_als)/w_aux;

return Rcpp::List::create(Rcpp::Named("w_aux") = w_aux,
                          Rcpp::Named("gamma") = gamma);

}
































































