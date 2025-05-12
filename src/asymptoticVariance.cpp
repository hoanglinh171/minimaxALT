#include "asymptoticVariance.h"


// Custom GSL error handler
void gsl_custom_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
    // std::cout << "GSL Error: " << reason << " at " << file << ":" << line << " (Error code: " << gsl_errno << ")" << std::endl;
    throw std::runtime_error(reason);
}


// Call this function before using GSL functions
void set_gsl_custom_handler() {
    gsl_set_error_handler(&gsl_custom_error_handler);
}


arma::vec softmax (arma::vec pi_arr) {
    int n_support = pi_arr.n_elem + 1;
    arma::vec pi_trans(n_support);

    double temp = 0;
    for (int i = 0; i < n_support - 1; i++) {
        pi_trans(i) = (1 - temp) * pi_arr(i);
        temp += pi_trans(i);
    }

    pi_trans(n_support - 1) = 1 - temp;

    return pi_trans;
}


arma::mat x_matrix(const arma::vec& alloc, int n_support, int n_factor, int n_unit) {
    if ((double)(alloc.n_elem + 1) / (double)n_support != (n_factor + 1)) {
        std::cout << "Number of supports and given allocation do not match." << std::endl;
        return arma::mat(1, 1, arma::fill::zeros);
    }

    arma::vec levels, pi_arr, level_pi;
    if (n_support <= 1) {
        levels = alloc;
        level_pi = {1};
    } else {
        levels = alloc.subvec(0, alloc.n_elem - n_support);
        pi_arr = alloc.subvec(alloc.n_elem - n_support + 1, alloc.n_elem - 1);
        level_pi = softmax(pi_arr);
    }

    // Calculate the number of units for each level
    level_pi = arma::round(level_pi * n_unit);
    level_pi(n_support - 1) = n_unit - sum(level_pi.head(n_support - 1));

    // Create the measure matrix with an intercept
    arma::mat measure_matrix(n_support, n_factor + 2);
    measure_matrix.col(0).ones();
    for (int i = 0; i < n_factor; ++i) {
        measure_matrix.col(i + 1) = levels.subvec(i * n_support, (i + 1) * n_support - 1);
    }
    measure_matrix.col(n_factor + 1) = level_pi;

    return measure_matrix;
}


arma::vec std_failure_time(double failure_time, const arma::mat& x_matrix, double sigma, const arma::vec& coef) {
    int n_var = coef.n_elem;
    int x_var = x_matrix.n_cols;
    int n_support = x_matrix.n_rows;
    if (n_var != x_var) {
        return {0};
    }

    arma::vec std_time(n_support, arma::fill::zeros);
    for (int i = 0; i < n_support; ++i) {
        double mu = arma::dot(x_matrix.row(i), coef);
        std_time(i) = (log(failure_time) - mu) / sigma;
    }
    return std_time;
}


arma::vec p_failure_time(const arma::vec& std_time, int distribution) {
    int n_support = std_time.n_elem;
    arma::vec phi_z(n_support, arma::fill::zeros);
    for (int i = 0; i < n_support; ++i) {
        if (distribution == 1) {
            phi_z(i) = 1 - exp(-exp(std_time(i))); // CDF of SEV distribution
        } else if (distribution == 2) {
            phi_z(i) = gsl_cdf_gaussian_P(std_time(i), 1); // CDF of normal distribution
        }
    }
    return phi_z;
}


double z_sqr_exp_z (double z, void * params) {
    (void)(params);
    return pow(log(z), 2) * z * exp(-z);
}



double z_exp_z (double z, void * params) {
    (void)(params);
    return log(z) * z * exp(-z);
}



double integrate_z_sqr_exp_z (double lower,
                             double upper) {

    gsl_set_error_handler(&gsl_custom_error_handler);

    double result, error;

    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(constants::INTEG_WORK_ALLOC);
    gsl_function F;
    F.function = &z_sqr_exp_z;
    F.params = nullptr;

    int status = gsl_integration_qags(&F, lower, upper, constants::INTEG_REL_TOL, constants::INTEG_ABS_TOL, constants::INTEG_WORK_ALLOC, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);

    if (status != GSL_SUCCESS) {
        throw std::runtime_error("GSL Integration Failed");
    }

    return result;
}



double integrate_z_exp_z (double lower,
                         double upper) {
    
    gsl_set_error_handler(&gsl_custom_error_handler);

    double result, error;

    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(constants::INTEG_WORK_ALLOC);
    gsl_function F;
    F.function = &z_exp_z;
    F.params = nullptr;

    int status = gsl_integration_qags(&F, lower, upper, constants::INTEG_REL_TOL, constants::INTEG_ABS_TOL, constants::INTEG_WORK_ALLOC, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);

    if (status != GSL_SUCCESS) {
        throw std::runtime_error("GSL Integration Failed");
    }

    return result;
}



double expected_hes_sigma(const arma::vec& std_censor,
                          const arma::vec& p_censor,
                          const arma::vec& alloc_unit,
                          double sigma,
                          int distribution) {
    int n_support = std_censor.n_elem;
    double n_info = 0;

    for (int i = 0; i < n_support; i++) {
        if (distribution == 1) {
            double expected_z_sqr_exp_z = integrate_z_sqr_exp_z(0, exp(std_censor[i]));
            double info_i = alloc_unit[i]
                            * (1 / pow(sigma, 2))
                            * (p_censor[i] + expected_z_sqr_exp_z
                               + (1 - p_censor[i]) * pow(std_censor[i], 2) * exp(std_censor[i]));
            n_info += info_i;
        } else if (distribution == 2) {
            double pdf_censor = gsl_ran_gaussian_pdf(std_censor[i], 1);
            double info_i = alloc_unit[i]
                            * (1 / pow(sigma, 2))
                            * (2 * p_censor[i]
                               - std_censor[i] * pdf_censor * (1 + pow(std_censor[i], 2) - std_censor[i] * pdf_censor / (1 - p_censor[i])));
            n_info += info_i;
        }
    }

    return n_info;
}



double expected_hes_beta_jk(int j,
                            int k,
                            const arma::mat& X,
                            const arma::vec& std_censor,
                            const arma::vec& p_censor,
                            const arma::vec& alloc_unit,
                            double sigma,
                            int distribution) {
    int n_support = std_censor.n_elem;
    double n_info = 0;

    for (int i = 0; i < n_support; i++) {
        if (distribution == 1) {
            double info_i = alloc_unit(i) * (X(i, j) * X(i, k)) / pow(sigma, 2) * p_censor(i);
            n_info += info_i;
        } else if (distribution == 2) {
            double pdf_censor = gsl_ran_gaussian_pdf(std_censor(i), 1);
            double info_i = alloc_unit(i)
                            * (X(i, j) * X(i, k)) / pow(sigma, 2)
                            * (p_censor(i) - pdf_censor * (std_censor(i) - pdf_censor / (1 - p_censor(i))));
            n_info += info_i;
        }
    }

    return n_info;
}



double expected_hes_sigma_beta_j(int j,
                                 const arma::mat& X,
                                 const arma::vec& std_censor,
                                 const arma::vec& p_censor,
                                 const arma::vec& alloc_unit,
                                 double sigma,
                                 int distribution) {
    int n_support = std_censor.n_elem;
    double n_info = 0;

    for (int i = 0; i < n_support; i++) {
        if (distribution == 1) {
            double expected_z_exp_z = integrate_z_exp_z(0, exp(std_censor(i)));
            double info_i = alloc_unit(i)
                            * X(i, j) / pow(sigma, 2)
                            * ((1 - p_censor[i]) * std_censor[i] * exp(std_censor[i]) + expected_z_exp_z);
            n_info += info_i;
        } else if (distribution == 2) {
            double pdf_censor = gsl_ran_gaussian_pdf(std_censor[i], 1);
            double info_i = alloc_unit[i]
                            * X(i, j) / pow(sigma, 2)
                            * (-pdf_censor * (1 + std_censor[i] * (std_censor[i] - pdf_censor / (1 - p_censor[i]))));
            n_info += info_i;
        }
    }

    return n_info;
}


arma::mat fisher_info(const arma::mat& X,
                      const arma::vec& std_censor,
                      const arma::vec& p_censor,
                      const arma::vec& alloc_unit,
                      double sigma,
                      int distribution) {
    int Nparams = X.n_cols + 1;

    arma::mat fisher_info_mat(Nparams, Nparams, arma::fill::zeros);

    for (int i = 0; i < Nparams; i++) {
        if (i == 0) {
            double hes_sigma = expected_hes_sigma(std_censor, p_censor, alloc_unit, sigma, distribution);
            fisher_info_mat(i, i) = hes_sigma;
        } else {
            double hes_beta_ii = expected_hes_beta_jk(i - 1, i - 1, X, std_censor, p_censor, alloc_unit, sigma, distribution);
            fisher_info_mat(i, i) = hes_beta_ii;
        }

        for (int j = i + 1; j < Nparams; j++) {
            if (i == 0) {
                double hes_sigma_betaj = expected_hes_sigma_beta_j(j - 1, X, std_censor, p_censor, alloc_unit, sigma, distribution);
                fisher_info_mat(i, j) = hes_sigma_betaj;
                fisher_info_mat(j, i) = hes_sigma_betaj;
            } else {
                double hes_beta_ij = expected_hes_beta_jk(i - 1, j - 1, X, std_censor, p_censor, alloc_unit, sigma, distribution);
                fisher_info_mat(i, j) = hes_beta_ij;
                fisher_info_mat(j, i) = hes_beta_ij;
            }
        }
    }

    return fisher_info_mat;
}


double is_singular(const arma::mat& matrix) {
    double cond_num = arma::rcond(matrix);
    if (cond_num <= constants::COND_NUM_LIM) {
        // std::cout << cond_num << std::endl;
        return true;
    }
    return false;
}



arma::mat cov_matrix(const arma::mat& fisher_info_mat) {
    int Nparams = fisher_info_mat.n_rows;
    arma::mat cov_matrix_mat(Nparams, Nparams, arma::fill::zeros);

    if (is_singular(fisher_info_mat)) {
        return cov_matrix_mat;
    } else {
        cov_matrix_mat = arma::inv(fisher_info_mat);
        return cov_matrix_mat;
    }
}


double inverse_cdf(double p,
                   int distribution) {
    double inv_cdf = 0;
    if (distribution == 1) {
        if (p >= constants::UPP_P) {
            // a - b = 3.47
            inv_cdf = log(-log(constants::LOW_P));
        } else if (p <= constants::LOW_P){
            inv_cdf = log(-log(constants::UPP_P));
        } else {
            inv_cdf = log(-log(1-p));
        }
    } else if (distribution == 2) {  // inv_cdf appx 7.
        if (p >=  constants::UPP_P) {
            p = constants::UPP_P;
        } else if (p <= constants::LOW_P) {
            p = constants::LOW_P;
        }
        inv_cdf = gsl_cdf_gaussian_Pinv(p, 1);
    }
    return inv_cdf;
}


double cdf(double z,
            int distribution) {
    double p = 0;
    if (distribution == 1) {
        p = 1 - exp(-exp(z));
    } else if (distribution == 2) {
        p = gsl_cdf_gaussian_P(z, 1);
    }
    return p;
}


double log_tp (const arma::vec& theta,
              double p,
              const arma::vec& use_cond,
              int distribution) {
    int n_theta = theta.n_elem;
    double log_tp = 0;
    log_tp = inverse_cdf(p, distribution) * theta(0);
    log_tp += arma::dot(theta.subvec(1, n_theta - 1), use_cond);

    return log_tp;
}



double gsl_log_tp(double theta_i, void* params) {
    deriv_params* dparams = static_cast<deriv_params*>(params);
    arma::uword theta_idx = dparams->theta_idx;
    arma::vec theta = dparams->theta;

    double original_value = theta(theta_idx);
    theta(theta_idx) = theta_i;
    double result = log_tp(theta, dparams->p, dparams->std_use_cond, dparams->distribution);
    theta(theta_idx) = original_value;

    return result;
}


arma::vec gradient(const arma::vec& theta,
                   double p,
                   const arma::vec& use_cond,
                   int distribution) {

    gsl_set_error_handler(&gsl_custom_error_handler);

    arma::uword n_theta = theta.n_elem;
    gsl_function F;
    arma::vec theta_vect = theta;
    double result, abserr;

    deriv_params params = {theta_vect, 0, p, use_cond, distribution};
    arma::vec deriv (n_theta);

    for (arma::uword i = 0; i < n_theta; ++i) {
        params.theta_idx = i;
        F.function = &gsl_log_tp;
        F.params = &params;
        gsl_deriv_central(&F, theta_vect(i), constants::DERIV_TOL, &result, &abserr);
        deriv(i) = result;
    }

    return deriv;
}


arma::vec reparameterize(double censor_time, const arma::vec p_coef, double sigma, int distribution, double x_l) {
    arma::uword Ncoef = p_coef.n_elem;
    arma::vec coef = p_coef;
    coef(0) = log(censor_time) - sigma * inverse_cdf(p_coef(0), distribution);
    for (arma::uword i = 1; i < Ncoef; ++i) {
        double inv_p1 = inverse_cdf(p_coef(i-1), distribution);
        double inv_p2 = inverse_cdf(p_coef(i), distribution);
        coef(i) = sigma * (inv_p1 - inv_p2) / (1 - x_l);
        coef(0) = coef(0) - coef(i) * x_l;
    }
    return coef;
}


arma::vec inverse_reparameterize(double censor_time, const arma::vec coef, double sigma, int distribution, double x_l) {
    arma::uword Ncoef = coef.n_elem;
    arma::vec p_coef = coef;
    p_coef(0) = (log(censor_time) - (coef(0) + x_l * arma::accu(coef.subvec(1, coef.n_elem - 1)))) / sigma;
    for (arma::uword i = 1; i < Ncoef; ++i) {
        p_coef(i) = coef(i) * (1 - x_l) / sigma;
        p_coef(i) = p_coef(i - 1) - p_coef(i);
    }

    for (arma::uword i = 0; i < Ncoef; ++i) {
        p_coef(i) = cdf(p_coef(i), distribution);
    }

    return p_coef;
}


int check_degenerate(arma::vec& p_coef,
                     design_info &design_info) {
    if (design_info.degenerate) {
        // Change number of factor
        design_info.n_factor = 1;

        // Change number of supports
        design_info.n_support = 2;

        // Change coefficient
        arma::vec coef(2);
        if (design_info.reparam) {
            coef = {p_coef(0),
                p_coef(p_coef.n_elem - 1)
            };
        } else {
            coef = {p_coef(0),
                arma::sum(p_coef.subvec(1, p_coef.n_elem - 1))
            };
        }
        p_coef.set_size(2);
        p_coef = coef;

        // Check use condition
        if (design_info.use_cond.n_elem > 1) {
            for (arma::uword i = 0; i < design_info.use_cond.n_elem; i++) {
                for (arma::uword j = i+1; j < design_info.use_cond.n_elem; j++) {
                    if (design_info.use_cond(i) != design_info.use_cond(j)) {
                        std::cout << "Degenerate design requires equal x's" << std::endl;
                        break;
                        return -1;
                    }
                }
            }
            design_info.use_cond.resize(1);
        }

    }
    return 0;
}

double obj_func(const arma::vec& alloc,
                const arma::vec& p_coef,
                int distribution,
                const design_info &design_info) {

    int n_support, n_factor, n_unit;
    arma::mat X, fisher_info_mat;
    double censor_time, sigma, p, x_l;
    arma::vec use_cond, coef, alloc_unit, std_censor, p_censor;
    std::string opt_type;

    n_support = design_info.n_support;
    n_factor = design_info.n_factor;
    n_unit = design_info.n_unit;
    censor_time = design_info.censor_time;
    sigma = design_info.sigma;
    p = design_info.p;
    use_cond = design_info.use_cond;
    opt_type = design_info.opt_type;
    x_l = design_info.x_l;
    coef = p_coef;


    // Reparameterize if needed
    if (design_info.reparam) {
        coef = reparameterize(censor_time, p_coef, sigma, distribution, x_l);
    }

    // Create measure matrix
    X = x_matrix(alloc, n_support, n_factor, n_unit);
    if (X.n_cols == 1) {
        return constants::BIG;
    }
    alloc_unit = X.col(X.n_cols - 1);
    X.shed_col(X.n_cols - 1);

    // Calculate standard censor time
    std_censor = std_failure_time(censor_time, X, sigma, coef);
    if (std_censor.n_elem == 1) {
        return constants::BIG;
    }

    // Calculate Phi(censor time)
    p_censor = p_failure_time(std_censor, distribution);

    // Calculate Fisher information matrix
    fisher_info_mat = fisher_info(X, std_censor, p_censor, alloc_unit, sigma, distribution);
    
    // Calculate covariance matrix
    arma::mat cov_mat = cov_matrix(fisher_info_mat);
    arma::vec variances = cov_mat.diag();
   

    if (is_singular(fisher_info_mat)) {
        return constants::BIG;

    }  else if (any(variances < 0)) {
        return constants::BIG;

    } else {
        arma::vec theta(coef.n_elem + 1);
        theta(0) = sigma;
        theta.subvec(1, theta.n_elem - 1) = coef;

        if (opt_type == "C") {  // C-optimality
            arma::mat use_cond_mtx = x_matrix(use_cond, 1, n_factor, 1);
            arma::vec use_cond_vec = arma::conv_to<arma::vec>::from(use_cond_mtx.row(0).subvec(0, theta.n_elem - 2));
            arma::vec deriv_log_tp = gradient(theta, p, use_cond_vec, distribution);

            // Calculate objective value
            double obj_value = arma::as_scalar(deriv_log_tp.t() * cov_mat * deriv_log_tp);
            return obj_value;

        } else if (opt_type == "D") {
            double det_fisher = arma::det(fisher_info_mat);
            return 1 / det_fisher;
        } else {
            std::cout << "Optimality type should be C or D." << std::endl;
            return constants::BIG;
        }
    }
}


double opt_crit(const arma::vec& alloc,
                const arma::vec& p_coef,
                int distribution,
                const design_info& design_info_in) {
    arma::vec coef = p_coef;
    design_info design_info_out = design_info_in;
    double obj_value;
    int check_degenerate_error;


    // Degenerate design requires equal x's
    check_degenerate_error = check_degenerate(coef, design_info_out);
    if (check_degenerate_error == -1) {
        return constants::BIG;
    }

    // Check condition of coefficient if reparam is true
    if (!design_info_out.reparam) {
        for (arma::uword i = 0; i < coef.n_elem; i++) {
            if (!std::isfinite(coef(i)) || coef(i) == 0) {
                return constants::BIG;
            }
        }

        coef = inverse_reparameterize(design_info_out.censor_time, coef, design_info_out.sigma, distribution, design_info_out.x_l);
    }

    // 0 <= p <= 1
    for (arma::uword i = 0; i < coef.n_elem; i++) {
        if (coef(i) < 0 || coef(i) > 1) {
            return constants::BIG;
        }
    }

    // pL <= pH
    for (arma::uword i = 0; i < coef.n_elem - 1; i++) {
        arma::uword j = i + 1;
        if (coef(i) >= coef(j)) {
            return constants::BIG;
        }
    }

    try {
        obj_value = obj_func(alloc, coef, distribution, design_info_out);
    } catch (const std::exception &e) {
        return constants::BIG;
    }

    if (!std::isfinite(obj_value)) {
        return constants::BIG;
    }

    return obj_value;
}
