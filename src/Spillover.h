#ifndef __Spillover__
#define __Spillover__

arma::vec rcpp_pgdraw(double b, 
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
                             double g_old);

double g_update(arma::mat x,
                arma::vec spillover_covar,
                arma::vec beta,
                double lambda,
                double alpha_g,
                double beta_g);

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
                      double acctot_phi_trans);

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
                        double acctot_theta_trans);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::vec spillover_covar,
                              arma::mat z, 
                              arma::vec beta,
                              double lambda,
                              arma::vec w);

Rcpp::List Spillover(int mcmc_samples,
                     arma::vec y,
                     arma::mat x,
                     arma::vec distance_to_ps,
                     arma::mat z,
                     arma::mat spatial_dists,
                     double a_theta_prior,
                     double b_theta_prior,
                     double metrop_var_phi_trans,
                     double metrop_var_theta_trans,
                     Rcpp::Nullable<double> alpha_g_prior,
                     Rcpp::Nullable<double> beta_g_prior,
                     Rcpp::Nullable<double> alpha_phi_prior,
                     Rcpp::Nullable<double> beta_phi_prior,
                     Rcpp::Nullable<double> alpha_sigma2_w_prior,
                     Rcpp::Nullable<double> beta_sigma2_w_prior,
                     Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                     Rcpp::Nullable<double> lambda_init,
                     Rcpp::Nullable<Rcpp::NumericVector> w_init,
                     Rcpp::Nullable<double> g_init,
                     Rcpp::Nullable<double> phi_init,
                     Rcpp::Nullable<double> theta_init,
                     Rcpp::Nullable<double> sigma2_w_init); 

#endif // __Spillover__
