#include "minimaxALT.h"

// [[Rcpp::export]]
Rcpp::List minimax_alt(int design_type, Rcpp::List &pso_info,
                       Rcpp::List &design_info_list,
                       Rcpp::List &init_bound_info,
                       Rcpp::List &nelder_mead_settings,
                       double n_threads, bool verbose) {
    
    #ifdef _OPENMP
    omp_set_num_threads(n_threads);
    #endif

    pso_options pso_opts;
    design_info design_info_local, design_info_glob;
    inner_optimization inner_param;
    pso_result pso_result_str;
    arma::mat init_coef_mat;

    get_pso_opts(pso_info, pso_opts);
    get_design_info(design_info_list, design_info_local, design_info_glob);
    get_inner_param(init_bound_info, inner_param);

    init_coef_mat = get_nelder_mead_settings(nelder_mead_settings);
    int n_multi_start = init_coef_mat.n_cols;

    if (verbose) {
        Rcpp::Rcout << "Calling Cpp PSO Kernel... " << std::endl;
    }

    pso_main(design_type, pso_opts, inner_param, design_info_local, design_info_glob, pso_result_str, n_multi_start, init_coef_mat, verbose);

    if (verbose) {
        Rcpp::Rcout << "Done." << std::endl;
    }

    arma::cube g_hist(pso_result_str.g_hist.n_cols, pso_result_str.g_hist.n_rows, pso_result_str.g_hist.n_slices);

    for (size_t i = 0; i < pso_result_str.g_hist.n_slices; ++i) {
        g_hist.slice(i) = pso_result_str.g_hist.slice(i).t();
    }

    return Rcpp::List::create(Rcpp::Named("g_best") = Rcpp::wrap(pso_result_str.g_best.t()),
                              Rcpp::Named("coef_best") = Rcpp::wrap(pso_result_str.coef_best.t()),
                              Rcpp::Named("distribution_best") = pso_result_str.distribution_best,
                              Rcpp::Named("max_directional_derivative") = pso_result_str.max_dd,
                              Rcpp::Named("fg_best") = pso_result_str.fg_best,
                              Rcpp::Named("fg_best_hist") = Rcpp::wrap(pso_result_str.fg_best_hist.t()),
                              Rcpp::Named("p_best") = Rcpp::wrap(pso_result_str.p_best.t()),
                              Rcpp::Named("fp_best") = Rcpp::wrap(pso_result_str.fp_best.t()),
                              Rcpp::Named("g_hist") = Rcpp::wrap(g_hist),
                              Rcpp::Named("coef_best_hist") = Rcpp::wrap(pso_result_str.coef_best_hist.t()),
                              Rcpp::Named("distribution_best_hist") = Rcpp::wrap(pso_result_str.distribution_best_hist.t()),
                              Rcpp::Named("model_set") = Rcpp::wrap(pso_result_str.model_set),
                              Rcpp::Named("model_weight") = Rcpp::wrap(pso_result_str.model_weight),
                              Rcpp::Named("equivalence_data") = Rcpp::wrap(pso_result_str.equivalence_data)
                        );
}


// [[Rcpp::export]]
Rcpp::NumericVector transform_proportion(Rcpp::NumericVector &dirichlet_prop) {
    arma::vec dirichlet_prop_temp(dirichlet_prop.begin(), dirichlet_prop.size(), false);
    arma::vec prop = softmax(dirichlet_prop_temp);

    return Rcpp::wrap(prop.t());
}


// [[Rcpp::export]]
Rcpp::NumericVector transform_sigmoid(Rcpp::NumericVector &inbound_sigmoid, 
                                    Rcpp::NumericVector &lower_bound, 
                                    Rcpp::NumericVector &upper_bound) {

    arma::vec inbound_sigmoid_temp(inbound_sigmoid.begin(), inbound_sigmoid.size(), false);
    arma::vec lower_bound_temp(lower_bound.begin(), lower_bound.size(), false);
    arma::vec upper_bound_temp(upper_bound.begin(), upper_bound.size(), false);

    arma::vec outbound_sigmoid = from_sigmoid(inbound_sigmoid_temp, lower_bound_temp, upper_bound_temp);

    return Rcpp::wrap(outbound_sigmoid.t());
}


// [[Rcpp::export]]
Rcpp::List equivalence_theorem(Rcpp::NumericVector &alloc, 
                                Rcpp::List &design_info_list, 
                                Rcpp::NumericMatrix &model_set_in) {
    
    pso_result pso_result_str;
    design_info design_info_glob, design_info_local;

    get_design_info(design_info_list, design_info_local, design_info_glob);
    
    arma::mat model_set (model_set_in.begin(), model_set_in.nrow(), model_set_in.ncol(), false);
    pso_result_str.model_set = model_set;

    arma::vec opt_alloc(alloc.begin(), alloc.size(), false);

    equivalence_plot_data(opt_alloc, design_info_glob, pso_result_str);

    return Rcpp::List::create(Rcpp::Named("max_directional_derivative") = pso_result_str.max_dd,
                              Rcpp::Named("model_set") = Rcpp::wrap(pso_result_str.model_set),
                              Rcpp::Named("model_weight") = Rcpp::wrap(pso_result_str.model_weight),
                              Rcpp::Named("equivalence_data") = Rcpp::wrap(pso_result_str.equivalence_data)
                        );

}
