#ifndef EQUIVALENCETHEOREM_H
#define EQUIVALENCETHEOREM_H

#include "common.h"
#include "asymptoticVariance.h"
#include "NelderMead.h"


arma::mat one_level_fisher(const arma::vec& opt_alloc, const arma::vec& p_coef,
                           const design_info &design_info_glob, arma::vec one_alloc_level);

double direction_deriv(const arma::vec& opt_alloc, const arma::vec& p_coef,
                       const design_info &design_info_glob, arma::vec one_alloc_level);

double sum_direction_derive(const arma::vec& gamma_weight, const arma::mat& gamma_set,
                            const arma::vec& opt_alloc, const design_info &design_info_glob, arma::vec one_alloc_level);

double lse_sum_der(const arma::vec& trans_gamma_weight, const arma::mat& gamma_set,
                   const arma::vec& opt_alloc, const design_info &design_info_glob,
                   double pct_lower, double pct_upper);

void lse_opt(const arma::vec &init_weigth, const arma::vec &opt_alloc,
             const design_info &design_info_glob, pso_result &pso_result_str);

arma::mat sim_sum_der(const arma::vec& gamma_weight, const arma::mat& gamma_set,
                      const arma::vec& opt_alloc, const design_info &design_info_glob);

void equivalence_plot_data(const arma::vec &opt_alloc,
                           const design_info &design_info_glob, pso_result &pso_result_str);


#endif // EQUIVALENCETHEOREM_H
