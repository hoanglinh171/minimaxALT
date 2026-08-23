#ifndef ASYMPTOTICVARIANCE_H
#define ASYMPTOTICVARIANCE_H

#include "common.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>


typedef struct {
    arma::vec theta;
    arma::uword theta_idx;
    double p;
    arma::vec std_use_cond;
    int distribution;
} deriv_params;


arma::vec softmax (arma::vec pi_arr);

arma::mat x_matrix(const arma::vec& alloc, int n_support, int n_factor, int n_unit);

arma::vec std_failure_time(double failure_time, const arma::mat& x_matrix, double sigma, const arma::vec& coef);

arma::vec p_failure_time(const arma::vec& std_time, int distribution);

double z_sqr_exp_z (double z, void * params);

double z_exp_z (double z, void * params);

double integrate_z_sqr_exp_z (double lower, double upper);

double integrate_z_exp_z (double lower, double upper);

double expected_hes_sigma(const arma::vec& std_censor, const arma::vec& p_censor, const arma::vec& alloc_unit,
                          double sigma, int distribution);
double expected_hes_beta_jk(int j,
                            int k,
                            const arma::mat& X,
                            const arma::vec& std_censor,
                            const arma::vec& p_censor,
                            const arma::vec& alloc_unit,
                            double sigma,
                            int distribution);

double expected_hes_sigma_beta_j(int j,
                                 const arma::mat& X,
                                 const arma::vec& std_censor,
                                 const arma::vec& p_censor,
                                 const arma::vec& alloc_unit,
                                 double sigma,
                                 int distribution);

arma::mat fisher_info(const arma::mat& X,
                      const arma::vec& std_censor,
                      const arma::vec& p_censor,
                      const arma::vec& alloc_unit,
                      double sigma,
                      int distribution);

double is_singular(const arma::mat& matrix);

arma::mat cov_matrix(const arma::mat& fisher_info_mat);

double log_tp (const arma::vec& theta,
              double p,
              const arma::vec& use_cond,
              int distribution);

double gsl_log_tp(double theta_i, void* params);

arma::vec gradient(const arma::vec& theta,
                   double p,
                   const arma::vec& use_cond,
                   int distribution);

double inverse_cdf(double p,
                   int distribution);

arma::vec reparameterize(double censor_time, const arma::vec p_coef, double sigma, int distribution, double x_l);

arma::vec inverse_reparameterize(double censor_time, const arma::vec coef, double sigma, int distribution, double x_l);

int check_degenerate(arma::vec& p_coef,
                     design_info &design_info);

double obj_func(const arma::vec& alloc,
                const arma::vec& p_coef,
                int distribution,
                const design_info &design_info);

double opt_crit(const arma::vec& alloc,
                const arma::vec& p_coef,
                int distribution,
                const design_info& design_info_in);

#endif // ASYMPTOTICVARIANCE_H
