#include "minimaxParams.h"

void get_pso_opts(Rcpp::List &pso_info, pso_options &pso_opts) {

    pso_opts.n_swarm = (int)Rcpp::as<int>(pso_info["n_swarm"]);
    pso_opts.d_swarm = (int)Rcpp::as<int>(pso_info["d_swarm"]);

    pso_opts.max_iter    = (int)Rcpp::as<int>(pso_info["max_iter"]);
    pso_opts.early_stopping    = (int)Rcpp::as<int>(pso_info["early_stopping"]);
    pso_opts.tol    = (double)Rcpp::as<double>(pso_info["tol"]);

    pso_opts.c1         = (double)Rcpp::as<double>(pso_info["c1"]);
    pso_opts.c2         = (double)Rcpp::as<double>(pso_info["c2"]);
    pso_opts.w0         = (double)Rcpp::as<double>(pso_info["w0"]);
    pso_opts.w1         = (double)Rcpp::as<double>(pso_info["w1"]);
    pso_opts.w_var      = (double)Rcpp::as<double>(pso_info["w_var"]);
    pso_opts.vk         = (double)Rcpp::as<double>(pso_info["vk"]);

    Rcpp::NumericVector var_upper_temp   = Rcpp::as<Rcpp::NumericVector>(pso_info["var_upper"]);
    arma::vec var_upper(var_upper_temp.begin(), var_upper_temp.size(), false);
    Rcpp::NumericVector var_lower_temp   = Rcpp::as<Rcpp::NumericVector>(pso_info["var_lower"]);
    arma::vec var_lower(var_lower_temp.begin(), var_lower_temp.size(), false);
    pso_opts.var_upper = var_upper;
    pso_opts.var_lower = var_lower;

    Rcpp::NumericMatrix init_swarm_temp = Rcpp::as<Rcpp::NumericMatrix>(pso_info["init_swarm"]);
    arma::mat init_swarm(init_swarm_temp.begin(), init_swarm_temp.nrow(), init_swarm_temp.ncol(), false);
    pso_opts.init_swarm  = init_swarm;

}


void get_design_info(Rcpp::List &design_info_list, design_info &design_info_local, design_info &design_info_glob) {

    // Minimax design info
    design_info_glob.n_support = (int)Rcpp::as<int>(design_info_list["n_support"]);
    design_info_glob.n_factor = (int)Rcpp::as<int>(design_info_list["n_factor"]);
    design_info_glob.n_unit = (int)Rcpp::as<int>(design_info_list["n_unit"]);

    design_info_glob.censor_time = (double)Rcpp::as<double>(design_info_list["censor_time"]);
    design_info_glob.sigma = (double)Rcpp::as<double>(design_info_list["sigma"]);
    design_info_glob.p = (double)Rcpp::as<double>(design_info_list["p"]);
    design_info_glob.x_l = (double)Rcpp::as<double>(design_info_list["x_l"]);
    design_info_glob.x_h = (double)Rcpp::as<double>(design_info_list["x_h"]);

    design_info_glob.reparam = (bool)Rcpp::as<bool>(design_info_list["reparam"]);
    design_info_glob.degenerate = (bool)Rcpp::as<bool>(design_info_list["degenerate"]);

    design_info_glob.opt_type = (std::string)Rcpp::as<std::string>(design_info_list["opt_type"]);

    Rcpp::NumericVector use_cond_temp  = Rcpp::as<Rcpp::NumericVector>(design_info_list["use_cond"]);
    arma::vec use_cond(use_cond_temp.begin(), use_cond_temp.size(), false);
    design_info_glob.use_cond = use_cond;


    // Local optimal design
    design_info_local.n_support = design_info_glob.n_support;
    design_info_local.n_factor = design_info_glob.n_factor;
    design_info_local.n_unit = design_info_glob.n_unit;

    design_info_local.censor_time = design_info_glob.censor_time;
    design_info_local.sigma = design_info_glob.sigma;
    design_info_local.p = design_info_glob.p;
    design_info_local.x_l = design_info_glob.x_l;
    design_info_local.x_h = design_info_glob.x_h;

    design_info_local.reparam = design_info_glob.reparam;
    design_info_local.degenerate = true;

    design_info_local.opt_type = design_info_glob.opt_type;

    design_info_local.use_cond = use_cond;
}


void get_inner_param(Rcpp::List &init_bound_info, inner_optimization &inner_param) {

    inner_param.model = (int)Rcpp::as<int>(init_bound_info["model"]);
    inner_param.opt_distribution = (int)Rcpp::as<int>(init_bound_info["opt_distribution"]);

    Rcpp::NumericVector coef_upper_temp   = Rcpp::as<Rcpp::NumericVector>(init_bound_info["coef_upper"]);
    arma::vec coef_upper(coef_upper_temp.begin(), coef_upper_temp.size(), false);
    Rcpp::NumericVector coef_lower_temp   = Rcpp::as<Rcpp::NumericVector>(init_bound_info["coef_lower"]);
    arma::vec coef_lower(coef_lower_temp.begin(), coef_lower_temp.size(), false);
    inner_param.coef_upper = coef_upper;
    inner_param.coef_lower = coef_lower;

    Rcpp::NumericVector init_coef_temp = Rcpp::as<Rcpp::NumericVector>(init_bound_info["init_coef"]);
    arma::vec init_coef(init_coef_temp.begin(), init_coef_temp.size(), false);
    inner_param.init_coef = init_coef;
    inner_param.opt_coef = to_sigmoid(init_coef, coef_lower, coef_upper);

    Rcpp::NumericVector local_upper_temp   = Rcpp::as<Rcpp::NumericVector>(init_bound_info["local_upper"]);
    arma::vec local_upper(local_upper_temp.begin(), local_upper_temp.size(), false);
    Rcpp::NumericVector local_lower_temp   = Rcpp::as<Rcpp::NumericVector>(init_bound_info["local_lower"]);
    arma::vec local_lower(local_lower_temp.begin(), local_lower_temp.size(), false);
    inner_param.local_upper = local_upper;
    inner_param.local_lower = local_lower;

    Rcpp::NumericVector init_local_temp = Rcpp::as<Rcpp::NumericVector>(init_bound_info["init_local"]);
    arma::vec init_local(init_local_temp.begin(), init_local_temp.size(), false);
    inner_param.init_local = init_local;
    inner_param.opt_local = to_sigmoid(init_local, local_lower, local_upper);

}


arma::mat get_nelder_mead_settings(Rcpp::List &nelder_mead_settings) {
    Rcpp::NumericMatrix init_coef_temp = Rcpp::as<Rcpp::NumericMatrix>(nelder_mead_settings["init_coef_mat"]);
    arma::mat init_coef_mat(init_coef_temp.begin(), init_coef_temp.nrow(), init_coef_temp.ncol(), false);
    
    return init_coef_mat;
}
