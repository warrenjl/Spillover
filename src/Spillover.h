#ifndef __Spillover__
#define __Spillover__

arma::vec rcpp_pgdraw(arma::vec b, 
                      arma::vec c);

Rcpp::List spatial_corr_fun(double phi,
                            arma::mat spatial_dists);

Rcpp::List w_aux_update(arma::vec y,
                        arma::mat x,
                        arma::vec spillover_covar,
                        arma::mat z,
                        arma::vec beta_old,
                        double lambda_old,
                        arma::vec w_old);

arma::vec beta_lambda_update(arma::mat x, 
                             arma::vec spillover_covar,
                             arma::mat z,
                             arma::vec w_aux,
                             arma::vec gamma,
                             arma::vec w_old,
                             double sigma2_regress);

arma::vec w_update(arma::mat x,
                   arma::vec spillover_covar,
                   arma::mat z,
                   arma::vec w_aux,
                   arma::vec gamma,
                   arma::vec beta,
                   double lambda,
                   double sigma2_w_old,
                   arma::mat corr_inv);

double sigma2_w_update(arma::vec w,
                       arma::mat corr_inv,
                       double alpha_sigma2_w,
                       double beta_sigma2_w);

Rcpp::List phi_update(arma::mat spatial_dists,
                      double phi_old,
                      arma::vec w,
                      double sigma2_w,
                      Rcpp::List spatial_corr_info,
                      double alpha_phi,
                      double beta_phi,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans);

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
                        int acctot_theta_trans);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::vec spillover_covar,
                              arma::mat z, 
                              arma::vec beta,
                              double lambda,
                              arma::vec w);

Rcpp::List Spillover(int mcmc_samples,
                     int spillover_covar_def, //1: Change point; 2: Exponential; 3: Gaussian 
                     arma::vec y,
                     arma::mat x,
                     arma::vec distance_to_ps,
                     arma::mat z,
                     arma::mat spatial_dists,
                     double metrop_var_phi_trans,
                     double metrop_var_theta_trans,
                     Rcpp::Nullable<double> sigma2_regress_prior,
                     Rcpp::Nullable<double> a_theta_prior,
                     Rcpp::Nullable<double> b_theta_prior,
                     Rcpp::Nullable<double> alpha_phi_prior,
                     Rcpp::Nullable<double> beta_phi_prior,
                     Rcpp::Nullable<double> alpha_sigma2_w_prior,
                     Rcpp::Nullable<double> beta_sigma2_w_prior,
                     Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                     Rcpp::Nullable<double> lambda_init,
                     Rcpp::Nullable<Rcpp::NumericVector> w_init,
                     Rcpp::Nullable<double> phi_init,
                     Rcpp::Nullable<double> theta_init,
                     Rcpp::Nullable<double> sigma2_w_init); 

#endif // __Spillover__
