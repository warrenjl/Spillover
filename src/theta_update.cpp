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
                        arma::vec spillover_covar_temp,
                        int spillover_covar_def,
                        double a_theta,
                        double b_theta,
                        double metrop_var_theta_trans,
                        int acctot_theta_trans){

int p_x = x.n_cols;
arma::vec beta_lambda(p_x + 1);
beta_lambda.subvec(0, (p_x - 1)) = beta;
beta_lambda(p_x) = lambda;
  
/*Second*/
double theta_trans_old = log((theta_old - a_theta)/(b_theta - theta_old));
arma::vec spillover_covar_old = (distance_to_ps <= theta_old)%spillover_covar_temp;
if(spillover_covar_def == 2){
  spillover_covar_old = (distance_to_ps <= theta_old)%exp(-distance_to_ps);
  }
if(spillover_covar_def == 3){
  spillover_covar_old = (distance_to_ps <= theta_old)%exp(-(distance_to_ps%distance_to_ps));
  }
arma::mat x_full_old = join_rows(x, spillover_covar_old);

double second = -0.50*dot((gamma - x_full_old*beta_lambda - z*w), (w_aux%(gamma - x_full_old*beta_lambda - z*w))) +
                theta_trans_old -
                2.00*log(1.00 + exp(theta_trans_old));

/*First*/
double theta_trans = R::rnorm(theta_trans_old, 
                              sqrt(metrop_var_theta_trans));
double theta = (b_theta*exp(theta_trans) + a_theta)/(1 + exp(theta_trans));
arma::vec spillover_covar = (distance_to_ps <= theta)%spillover_covar_temp;
if(spillover_covar_def == 2){
  spillover_covar = (distance_to_ps <= theta)%exp(-distance_to_ps);
  }
if(spillover_covar_def == 3){
  spillover_covar = (distance_to_ps <= theta)%exp(-(distance_to_ps%distance_to_ps));
  }
arma::mat x_full = join_rows(x, spillover_covar);

double first = -0.50*dot((gamma - x_full*beta_lambda - z*w), (w_aux%(gamma - x_full*beta_lambda - z*w))) +
               theta_trans -
               2.00*log(1.00 + exp(theta_trans));

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  theta = theta_old;
  acc = 0;
  }
acctot_theta_trans = acctot_theta_trans + 
                     acc;

return Rcpp::List::create(Rcpp::Named("theta") = theta,
                          Rcpp::Named("acctot_theta_trans") = acctot_theta_trans);

}
                 
  
