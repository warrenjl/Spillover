#include "RcppArmadillo.h"
#include "Spillover.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List theta_update(arma::mat x,
                        arma::mat z,
                        arma::vec distance_to_ps,
                        double theta_old,
                        arma::vec w_aux,
                        arma::vec gamma,
                        arma::vec beta,
                        double lambda,
                        arma::vec w,
                        double g,
                        double a_theta,
                        double b_theta,
                        double metrop_var_theta_trans,
                        double acctot_theta_trans){

arma::vec beta_lambda(x.n_cols + 1);
beta_lambda.subvec(0, (x.n_cols - 1)) = beta;
beta_lambda(x.n_cols) = lambda;
  
/*Second*/
double theta_trans_old = log((theta_old - a_theta)/(b_theta - theta_old));
arma::vec spillover_covar_old = (distance_to_ps <= theta_old)%exp(-(distance_to_ps%distance_to_ps));
arma::mat x_full_old = join_rows(x, spillover_covar_old);
double log_deter_old = 0; 
double sign_old = 0;     
log_det(log_deter_old, sign_old, inv_sympd(trans(x_full_old)*x_full_old));

double second = -0.50*dot((gamma - x_full_old*beta_lambda - z*w), (w_aux%(gamma - x_full_old*beta_lambda - z*w))) -
                0.50*log_deter_old -
                (0.50/g)*dot(beta_lambda, (trans(x_full_old)*x_full_old)*beta_lambda) +
                theta_trans_old -
                2*log(1 + exp(theta_trans_old));

/*First*/
double theta_trans = R::rnorm(theta_trans_old, 
                              sqrt(metrop_var_theta_trans));
double theta = (b_theta*exp(theta_trans) + a_theta)/(1 + exp(theta_trans));
arma::vec spillover_covar = (distance_to_ps <= theta)%exp(-(distance_to_ps%distance_to_ps));
arma::mat x_full = join_rows(x, spillover_covar);
double log_deter = 0; 
double sign = 0;     
log_det(log_deter, sign, inv_sympd(trans(x_full)*x_full));

double first = -0.50*dot((gamma - x_full*beta_lambda - z*w), (w_aux%(gamma - x_full*beta_lambda - z*w))) -
               0.50*log_deter -
               (0.50/g)*dot(beta_lambda, (trans(x_full)*x_full)*beta_lambda) + 
               theta_trans -
               2*log(1 + exp(theta_trans));

/*Decision*/
double ratio = exp(first - second);   
double acc = 1;
if(ratio < R::runif(0, 1)){
  theta = theta_old;
  acc = 0;
  }
acctot_theta_trans = acctot_theta_trans + 
                     acc;

return Rcpp::List::create(Rcpp::Named("theta") = theta,
                          Rcpp::Named("acctot_theta_trans") = acctot_theta_trans);

}
                 
  
