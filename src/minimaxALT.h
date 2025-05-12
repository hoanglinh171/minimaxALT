#ifndef MINIMAXALT_H
#define MINIMAXALT_H

#include "minimaxParams.h"
#include "asymptoticVariance.h"
#include "NelderMead.h"
#include "psoMain.h"
#include "equivalenceTheorem.h"

Rcpp::List minimax_alt(int design_type, Rcpp::List &pso_info,
    Rcpp::List &design_info_list,
    Rcpp::List &init_bound_info,
    Rcpp::List &nelder_mead_settings,
    double n_threads, bool verbose);

Rcpp::NumericVector transform_proportion(Rcpp::NumericVector &dirichlet_prop);

Rcpp::List equivalence_theorem(Rcpp::NumericVector &alloc, 
    Rcpp::List &design_info_list, 
    Rcpp::NumericMatrix &model_set_in);
    

#endif // MINIMAXALT_H
