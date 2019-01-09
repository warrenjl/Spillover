#include "RcppArmadillo.h"
#include "Spillover.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List Spillover(int mcmc_samples,
                     int spillover_covar_def, //1: Change point; 2: Exponential; 3: Gaussian 
                     arma::vec y,
                     arma::mat x,
                     arma::vec distance_to_ps,
                     arma::mat z,
                     arma::mat spatial_dists,
                     double metrop_var_phi_trans,
                     double metrop_var_theta_trans,
                     Rcpp::Nullable<double> sigma2_regress_prior = R_NilValue,
                     Rcpp::Nullable<double> a_theta_prior = R_NilValue,
                     Rcpp::Nullable<double> b_theta_prior = R_NilValue,
                     Rcpp::Nullable<double> alpha_phi_prior = R_NilValue,
                     Rcpp::Nullable<double> beta_phi_prior = R_NilValue,
                     Rcpp::Nullable<double> alpha_sigma2_w_prior = R_NilValue,
                     Rcpp::Nullable<double> beta_sigma2_w_prior = R_NilValue,
                     Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                     Rcpp::Nullable<double> lambda_init = R_NilValue,
                     Rcpp::Nullable<Rcpp::NumericVector> w_init = R_NilValue,
                     Rcpp::Nullable<double> phi_init = R_NilValue,
                     Rcpp::Nullable<double> theta_init = R_NilValue,
                     Rcpp::Nullable<double> sigma2_w_init = R_NilValue){
  
//Defining Parameters and Quantities of Interest
int n = y.size();
int p_x = x.n_cols;
int p_z = z.n_cols;
arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
arma::vec lambda(mcmc_samples); lambda.fill(0.00);
arma::mat w(p_z, mcmc_samples); w.fill(0.00);
arma::vec phi(mcmc_samples); phi.fill(0.00);
arma::vec theta(mcmc_samples); theta.fill(0.00);
arma::vec sigma2_w(mcmc_samples); sigma2_w.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);
  
//Prior Information
double sigma2_regress = 10000.00;
if(sigma2_regress_prior.isNotNull()){
  sigma2_regress = Rcpp::as<double>(sigma2_regress_prior);
  }

double a_theta = min(distance_to_ps);
if(a_theta_prior.isNotNull()){
  a_theta = Rcpp::as<double>(a_theta_prior);
  }

double b_theta = max(distance_to_ps);
if(b_theta_prior.isNotNull()){
  b_theta = Rcpp::as<double>(b_theta_prior);
  }
  
double alpha_phi = 1.00;
if(alpha_phi_prior.isNotNull()){
  alpha_phi = Rcpp::as<double>(alpha_phi_prior);
  }
  
double beta_phi = 1.00;
if(beta_phi_prior.isNotNull()){
  beta_phi = Rcpp::as<double>(beta_phi_prior);
  }

double alpha_sigma2_w = 3.00;
if(alpha_sigma2_w_prior.isNotNull()){
  alpha_sigma2_w = Rcpp::as<double>(alpha_sigma2_w_prior);
  }

double beta_sigma2_w = 2.00;
if(beta_sigma2_w_prior.isNotNull()){
  beta_sigma2_w = Rcpp::as<double>(beta_sigma2_w_prior);
  }

//Initial Values
beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

lambda(0) = 0.00;
if(lambda_init.isNotNull()){
  lambda(0) = Rcpp::as<double>(lambda_init);
  } 

w.col(0).fill(0.00);
if(w_init.isNotNull()){
  w.col(0) = Rcpp::as<arma::vec>(w_init);
  }

phi(0) = 0.01*max(distance_to_ps);
if(phi_init.isNotNull()){
  phi(0) = Rcpp::as<double>(phi_init);
  }

theta(0) = 0.50*(b_theta - a_theta);
if(theta_init.isNotNull()){
  theta(0) = Rcpp::as<double>(theta_init);
  }

sigma2_w(0) = 1.00;
if(sigma2_w_init.isNotNull()){
  sigma2_w(0) = Rcpp::as<double>(sigma2_w_init);
  }

Rcpp::List spatial_corr_info = spatial_corr_fun(phi(0),
                                                spatial_dists);

arma::vec spillover_covar_temp(n); spillover_covar_temp.fill(1.00);
arma::vec spillover_covar = (distance_to_ps <= theta(0))%spillover_covar_temp;
if(spillover_covar_def == 2){
  spillover_covar = (distance_to_ps <= theta(0))%exp(-distance_to_ps);
  }
if(spillover_covar_def == 3){
  spillover_covar = (distance_to_ps <= theta(0))%exp(-(distance_to_ps%distance_to_ps));
  }

neg_two_loglike(0) = neg_two_loglike_update(y,
                                            x,
                                            spillover_covar,
                                            z, 
                                            beta.col(0),
                                            lambda(0),
                                            w.col(0));

//Metropolis Settings
int acctot_phi_trans = 0;
int acctot_theta_trans = 0;

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++j){
  
   //w_aux Update
   Rcpp::List w_aux_output = w_aux_update(y,
                                          x,
                                          spillover_covar,
                                          z,
                                          beta.col(j-1),
                                          lambda(j-1),
                                          w.col(j-1));
  
   arma::vec w_aux = w_aux_output[0];
   arma::vec gamma = w_aux_output[1];
  
   //beta_lambda Update
   arma::vec beta_lambda = beta_lambda_update(x, 
                                              spillover_covar,
                                              z,
                                              w_aux,
                                              gamma,
                                              w.col(j-1),
                                              sigma2_regress);
   
   beta.col(j) = beta_lambda.subvec(0, (p_x - 1));
   lambda(j) = beta_lambda(p_x);
  
   //w Update
   w.col(j) = w_update(x,
                       spillover_covar,
                       z,
                       w_aux,
                       gamma,
                       beta.col(j),
                       lambda(j),
                       sigma2_w(j-1),
                       spatial_corr_info[0]);
   
   //sigma2_w update
   sigma2_w(j) = sigma2_w_update(w.col(j),
                                 spatial_corr_info[0],
                                 alpha_sigma2_w,
                                 beta_sigma2_w);
   
   //phi Update
   Rcpp::List phi_output = phi_update(spatial_dists,
                                      phi(j-1),
                                      w.col(j),
                                      sigma2_w(j),
                                      spatial_corr_info,
                                      alpha_phi,
                                      beta_phi,
                                      metrop_var_phi_trans,
                                      acctot_phi_trans);
     
   phi(j) = phi_output[0];
   acctot_phi_trans = phi_output[1];
   spatial_corr_info = phi_output[2];
   
   //theta Update
   Rcpp::List theta_output = theta_update(x,
                                          z,
                                          distance_to_ps,
                                          theta(j-1),
                                          w_aux,
                                          gamma,
                                          beta.col(j),
                                          lambda(j),
                                          w.col(j),
                                          spillover_covar_temp,
                                          spillover_covar_def,
                                          a_theta,
                                          b_theta,
                                          metrop_var_theta_trans,
                                          acctot_theta_trans);
   
   theta(j) = theta_output[0];
   acctot_theta_trans = theta_output[1];
   spillover_covar = (distance_to_ps <= theta(j))%spillover_covar_temp;
   if(spillover_covar_def == 2){
     spillover_covar = (distance_to_ps <= theta(j))%exp(-distance_to_ps);
     }
   if(spillover_covar_def == 3){
     spillover_covar = (distance_to_ps <= theta(j))%exp(-(distance_to_ps%distance_to_ps));
     }
   
   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update(y,
                                               x,
                                               spillover_covar,
                                               z, 
                                               beta.col(j),
                                               lambda(j),
                                               w.col(j));
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
     
   if(((j + 1) % int(round(mcmc_samples*0.10)) == 0)){
     if(spillover_covar_def == 1){
       Rcpp::Rcout << "**********************" << std::endl;
       Rcpp::Rcout << "Change Point Spillover" << std::endl;
       }
     if(spillover_covar_def == 2){
       Rcpp::Rcout << "*********************" << std::endl;
       Rcpp::Rcout << "Exponential Spillover" << std::endl;
       }
     if(spillover_covar_def == 3){
       Rcpp::Rcout << "*********************" << std::endl;
       Rcpp::Rcout << "Gaussian Spillover" << std::endl;
      }
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     double accrate_phi_trans = round(100*(acctot_phi_trans/(double)j));
     Rcpp::Rcout << "phi Acceptance: " << accrate_phi_trans << "%" << std::endl;
     double accrate_theta_trans = round(100*(acctot_theta_trans/(double)j));
     Rcpp::Rcout << "theta Acceptance: " << accrate_theta_trans << "%" << std::endl;
     }
     
   }
  
return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("lambda") = lambda,
                          Rcpp::Named("w") = w,
                          Rcpp::Named("sigma2_w") = sigma2_w,
                          Rcpp::Named("phi") = phi,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans,
                          Rcpp::Named("acctot_theta_trans") = acctot_theta_trans);

}
