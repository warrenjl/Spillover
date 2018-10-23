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
                        double a_theta,
                        double b_theta,
                        double metrop_var_theta_trans,
                        double acctot_theta_trans){

/*Second*/
double theta_trans_old = log((theta_old - a_theta)/(b_theta - theta_old));
arma::vec spillover_covar_old = (distance_to_ps <= theta_old)%exp(-(distance_to_ps%distance_to_ps));

double second = -0.50*dot((gamma - x*beta - spillover_covar_old*lambda - z*w), (w_aux%(gamma - x*beta - spillover_covar_old*lambda - z*w))) + 
                theta_trans_old -
                2*log(1 + exp(theta_trans_old));

/*First*/
double theta_trans = R::rnorm(theta_trans_old, 
                              sqrt(metrop_var_theta_trans));
double theta = (b_theta*exp(theta_trans) + a_theta)/(1 + exp(theta_trans));
arma::vec spillover_covar = (distance_to_ps <= theta)%exp(-(distance_to_ps%distance_to_ps));

double first = -0.50*dot((gamma - x*beta - spillover_covar*lambda - z*w), (w_aux%(gamma - x*beta - spillover_covar*lambda - z*w))) +  
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
                 
  
