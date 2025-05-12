#ifndef MINIMAXPARAMS_H
#define MINIMAXPARAMS_H

#include "common.h"

void get_pso_opts(Rcpp::List &pso_info, pso_options &pso_opts);
void get_design_info(Rcpp::List &design_info_list, design_info &design_info_local, design_info &design_info_glob);
void get_inner_param(Rcpp::List &init_bound_info, inner_optimization &inner_param);
arma::mat get_nelder_mead_settings(Rcpp::List &nelder_mead_settings);

#endif // MINIMAXPARAMS_H
