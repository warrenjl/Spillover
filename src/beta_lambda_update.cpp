#include "RcppArmadillo.h"
#include "Spillover.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec beta_lambda_update(arma::mat x, 
                             arma::vec spillover_covar,
                             arma::mat z,
                             arma::vec w_aux,
                             arma::vec gamma,
                             arma::vec w_old,
                             double sigma2_regress){

arma::mat x_full = join_rows(x, spillover_covar);
int p_x_full = x_full.n_cols;
int n = w_aux.size();

arma::mat w_aux_mat(n, p_x_full);
for(int j = 0; j < p_x_full; ++j){
   w_aux_mat.col(j) = w_aux;
   }

arma::mat x_full_trans = trans(x_full);

arma::mat cov_beta_lambda = inv_sympd(x_full_trans*(w_aux_mat%x_full) + 
                                      eye(p_x_full, p_x_full)/sigma2_regress);

arma::vec mean_beta_lambda = cov_beta_lambda*(x_full_trans*(w_aux%(gamma - z*w_old)));

arma::mat ind_norms = arma::randn(1, p_x_full);
arma::vec beta_lambda = mean_beta_lambda + 
                        trans(ind_norms*arma::chol(cov_beta_lambda));

return beta_lambda;

}



