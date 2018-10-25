#include "RcppArmadillo.h"
#include "Spillover.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double g_update(arma::mat x,
                arma::vec spillover_covar,
                arma::vec beta,
                double lambda,
                double alpha_g,
                double beta_g){

int p_x_full = x.n_cols + 1;

arma::vec beta_lambda(x.n_cols + 1);
beta_lambda.subvec(0, (x.n_cols - 1)) = beta;
beta_lambda(x.n_cols) = lambda;

arma::mat x_full = join_rows(x, spillover_covar);

double alpha_g_update = 0.50*p_x_full + 
                        alpha_g;

double beta_g_update = 0.50*dot(beta_lambda, ((trans(x_full)*x_full)*beta_lambda)) + 
                       beta_g;

double g = 1/R::rgamma(alpha_g_update,
                       (1/beta_g_update));

return(g);

}





