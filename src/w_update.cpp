#include "RcppArmadillo.h"
#include "Spillover.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec w_update(arma::mat x,
                   arma::vec spillover_covar,
                   arma::mat z,
                   arma::vec w_aux,
                   arma::vec gamma,
                   arma::vec beta,
                   double lambda,
                   double sigma2_w_old,
                   arma::mat corr_inv){

int p_z = z.n_cols;
int n = w_aux.size();

arma::mat w_aux_mat(n, p_z);
for(int j = 0; j < p_z; ++j){
   w_aux_mat.col(j) = w_aux;
   }

arma::mat z_trans = trans(z);
  
arma::mat cov_w = inv_sympd(z_trans*(w_aux_mat%z) + 
                            (1/sigma2_w_old)*corr_inv);

arma::vec mean_w = cov_w*(z_trans*(w_aux%(gamma - x*beta - spillover_covar*lambda)));
  
arma::mat ind_norms = arma::randn(1, p_z);
arma::vec w = mean_w + 
              trans(ind_norms*arma::chol(cov_w));

w = w - mean(w);
  
return(w);
  
}

