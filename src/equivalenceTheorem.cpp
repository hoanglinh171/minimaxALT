#include "equivalenceTheorem.h"
#include "NelderMead.h"


// Function to calculate one-level Fisher information matrix
arma::mat one_level_fisher(const arma::vec& opt_alloc, const arma::vec& p_coef, int distribution,
                           const design_info &design_info_glob, arma::vec one_alloc_level) {

    int n_support = design_info_glob.n_support;
    int n_factor = design_info_glob.n_factor;
    int n_unit = design_info_glob.n_unit;
    double censor_time = design_info_glob.censor_time;
    double sigma = design_info_glob.sigma;
    double x_l = design_info_glob.x_l;
    arma::vec coef;                        
    
    if (design_info_glob.reparam) {
        coef = reparameterize(censor_time, p_coef, sigma, distribution, x_l);
    } else {
        coef = p_coef;
    }

    arma::mat X = x_matrix(opt_alloc, n_support, n_factor, n_unit);
    X.shed_col(X.n_cols - 1);

    // Let all levels have 0 unit, except the level at one_alloc_idx
    arma::vec alloc_unit = arma::zeros(n_support);
    alloc_unit(0) = n_unit;
    X(0, arma::span(1, n_factor)) = one_alloc_level.t();

    arma::vec std_censor = std_failure_time(censor_time, X, sigma, coef);
    arma::vec p_censor = p_failure_time(std_censor, distribution);

    arma::mat fisher_mtx = fisher_info(X, std_censor, p_censor, alloc_unit, sigma, distribution);
    return fisher_mtx;
}


// Directional derivative of one-level design
double direction_deriv(const arma::vec& opt_alloc, const arma::vec& model_vec,
                       const design_info &design_info_glob, arma::vec one_alloc_level) {
    int n_support = design_info_glob.n_support;
    int n_factor = design_info_glob.n_factor;
    int n_unit = design_info_glob.n_unit;
    double censor_time = design_info_glob.censor_time;
    double sigma = design_info_glob.sigma;
    double x_l = design_info_glob.x_l;
    double p = design_info_glob.p;
    arma::vec use_cond = design_info_glob.use_cond;

    int distribution = model_vec(model_vec.size() - 1);
    arma::vec p_coef = model_vec.subvec(0, model_vec.size() - 2);

    arma::vec coef = reparameterize(censor_time, p_coef, sigma, distribution, x_l);
    arma::mat X = x_matrix(opt_alloc, n_support, n_factor, n_unit);
    arma::vec alloc_unit = X.col(X.n_cols - 1);
    X.shed_col(X.n_cols - 1);

    arma::vec std_censor = std_failure_time(censor_time, X, sigma, coef);
    arma::vec p_censor = p_failure_time(std_censor, distribution);

    arma::mat fisher_mtx = fisher_info(X, std_censor, p_censor, alloc_unit, sigma, distribution);
    arma::mat cov_mtx = cov_matrix(fisher_mtx);

    arma::mat fisher_mtx_one = one_level_fisher(opt_alloc, p_coef, distribution, design_info_glob, one_alloc_level);

    arma::vec deriv_log_tp = gradient(arma::join_cols(arma::vec({sigma}), coef), p,
                                      arma::join_cols(arma::vec({1}), use_cond), distribution);

    double nom = arma::as_scalar(deriv_log_tp.t() * cov_mtx * fisher_mtx_one * cov_mtx * deriv_log_tp);
    double denom = arma::as_scalar(deriv_log_tp.t() * cov_mtx * deriv_log_tp);

    return nom / denom;
}


double sum_direction_derive(const arma::vec& model_weight, const arma::mat& model_set,
                            const arma::vec& opt_alloc, const design_info &design_info_glob, arma::vec one_alloc_level) {
    double sum_deriv = 0.0;


    for (size_t i = 0; i < model_weight.n_elem; ++i) {
        // Extract the i-th row of model_set as a vector
        arma::rowvec model_set_row = model_set.row(i);
        arma::vec model_vec = model_set_row.t(); // Transpose row to column vector

        // Call the direction_deriv function
        double dir_deriv = direction_deriv(opt_alloc, model_vec, design_info_glob, one_alloc_level);

        // Accumulate the weighted directional derivative
        sum_deriv += model_weight(i) * dir_deriv;
    }

    return sum_deriv;
}


double lse_sum_der(const arma::vec& trans_model_weight, const arma::mat& model_set,
                   const arma::vec& opt_alloc, const design_info &design_info_glob,
                   double pct_lower, double pct_upper) {
    // Ensure sum of weights = 1 using softmax
    arma::vec weigth_lower(trans_model_weight.n_elem); weigth_lower.fill(pct_lower);
    arma::vec weigth_upper(trans_model_weight.n_elem); weigth_upper.fill(pct_upper);
    arma::vec model_weight = to_sigmoid(trans_model_weight, weigth_lower, weigth_upper);
    model_weight = softmax(model_weight);

    double n_support = design_info_glob.n_support;
    double n_factor = design_info_glob.n_factor;
    double sum_deriv;

    double obj_value = 0.0;
    arma::mat X = x_matrix(opt_alloc, n_support, n_factor,
                           design_info_glob.n_unit);
    arma::vec alloc_unit = X.col(X.n_cols - 1);

    // Loop through each support level
    for (int i = 0; i < n_support; ++i) {
        // Calculate the sum of directional derivatives
        if (alloc_unit(i) > 0) {
            arma::vec one_alloc_level = X(i, arma::span(1, n_factor)).t();
            sum_deriv = sum_direction_derive(model_weight, model_set, opt_alloc,
                                             design_info_glob, one_alloc_level);

            // Each derivative equals 1 if the design is minimax
            obj_value += std::pow(sum_deriv - 1.0, 2);
        }
    }

    return obj_value;
}


void lse_opt(const arma::vec &init_weigth, const arma::vec &opt_alloc,
                  const design_info &design_info_glob, pso_result &pso_result_str) {
    arma::mat model_set = pso_result_str.model_set;
    int n = model_set.n_rows - 1;

    arma::vec par, par_i;
    double sum_der;
    double upp = 1, low = 0;
    arma::vec weigth_lower(n); weigth_lower.fill(low);
    arma::vec weigth_upper(n); weigth_upper.fill(upp);

    if (model_set.n_rows == 1) {
        par = {1};
    } else {
        std::tie(par, sum_der) = nelder_mead(NM_PARAMS, init_weigth, lse_sum_der, model_set, opt_alloc, design_info_glob, low, upp);
        par = softmax(to_sigmoid(par, weigth_lower, weigth_upper));
    }
    
    pso_result_str.model_weight = par;

}


void generate_combinations(const arma::mat& vectors, arma::mat& combinations,
                           arma::rowvec& temp, int depth, int& index) {
    if (depth == (int)vectors.n_cols) {
        combinations.row(index) = temp;
        index++;
        return;
    }

    for (arma::uword i = 0; i < vectors.n_rows; ++i) {
        temp(depth) = vectors(i, depth);
        generate_combinations(vectors, combinations, temp, depth + 1, index);
    }
}


arma::mat sim_sum_der(const arma::vec& model_weight, const arma::mat& model_set,
                      const arma::vec& opt_alloc, const design_info &design_info_glob) {
    int n_factor = design_info_glob.n_factor;
    int sim_num = 101;
    int sim_elem_num = pow(sim_num, n_factor);
    arma::mat vectors(sim_num, n_factor);
    vectors.each_col() = arma::linspace(0, 1, sim_num);
    arma::mat combinations(sim_elem_num, n_factor);
    arma::rowvec temp(n_factor);
    int index = 0;

    generate_combinations(vectors, combinations, temp, 0, index);
    arma::vec dd_lst(sim_elem_num);

    for (int i = 0; i < sim_elem_num; ++i) {
        arma::vec one_alloc_level = combinations.row(i).t();
        double sum_dd = sum_direction_derive(model_weight, model_set, opt_alloc, design_info_glob, one_alloc_level);
        dd_lst(i) = sum_dd;
    }

    arma::mat data(sim_elem_num, n_factor + 2);
    data(arma::span::all, arma::span(0, data.n_cols - 3)) = combinations;
    data.col(data.n_cols - 2) = dd_lst;
    data.col(data.n_cols - 1).fill(0);

    return data;
}


void equivalence_plot_data(const arma::vec &opt_alloc,
                                const design_info &design_info_glob, pso_result &pso_result_str) {
    arma::mat model_set = pso_result_str.model_set, sim_data, plot_data, max_plot_data;
    int n = model_set.n_rows - 1;

    int n_support = design_info_glob.n_support;
    int n_factor = design_info_glob.n_factor;
    // double x_l = design_info_glob.x_l, x_h = design_info_glob.x_h;
    double opt_level_dd;
    arma::vec init_weigth(n), model_weigth(n);

    arma::mat level_data(n_support, n_factor + 2);
    level_data.col(level_data.n_cols - 1).fill(1);
    arma::mat X = x_matrix(opt_alloc, n_support, n_factor,
                           design_info_glob.n_unit);

    double max_dd = 10000, max_dd_i;
    arma::vec optimal_weigth;

    for (int i = 0; i < 3; i++) {
        init_weigth = 2 * arma::randu(n);
        lse_opt(init_weigth, opt_alloc, design_info_glob, pso_result_str);
        model_weigth = pso_result_str.model_weight;


        for (int i = 0; i < n_support; i++) {
            arma::vec one_alloc_level = X(i, arma::span(1, n_factor)).t();
            opt_level_dd = sum_direction_derive(model_weigth, model_set, opt_alloc, design_info_glob, one_alloc_level);
            level_data(i, arma::span(0, n_factor - 1)) = X(i, arma::span(1, n_factor));
            level_data(i, level_data.n_cols - 2) = opt_level_dd;
        }

        sim_data = sim_sum_der(model_weigth, model_set, opt_alloc, design_info_glob);
        plot_data = arma::join_cols(level_data, sim_data);
        
        // Extract first (n-1) columns
        arma::mat plot_data_sub = plot_data.cols(0, plot_data.n_cols - 3);
        
        // Logical mask: each row where all values in [x_l, x_h]
        arma::uvec valid_rows = arma::all((plot_data_sub >= design_info_glob.x_l) && (plot_data_sub <= design_info_glob.x_h), 1);  
        
        // Filter corresponding n-th column values
        arma::vec nth_col_valid = plot_data.col(plot_data.n_cols - 2);
        nth_col_valid = nth_col_valid.elem(arma::find(valid_rows));
        
        max_dd_i = nth_col_valid.max();
        if (abs(max_dd_i - 1) < abs(max_dd - 1)) {
            max_dd = max_dd_i;
            optimal_weigth = model_weigth;
            max_plot_data = plot_data;
        }

    }

    pso_result_str.model_weight = optimal_weigth;
    pso_result_str.max_dd = max_dd;
    pso_result_str.equivalence_data = max_plot_data;

}

