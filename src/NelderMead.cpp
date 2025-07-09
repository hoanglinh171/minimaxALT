#include "NelderMead.h"


double avar_wrapper(const arma::vec &init_local, inner_optimization &inner_param,
                    const design_info &design_info_local, int distribution) {
    arma::vec local_sigmoid(init_local.n_elem);
    double avar;

    local_sigmoid = to_sigmoid(init_local, inner_param.local_lower, inner_param.local_upper);

    avar = opt_crit(local_sigmoid, inner_param.init_coef, distribution, design_info_local);

    return avar;
}


double local_optimal(inner_optimization &inner_param, const design_info &design_info_local, int distribution) {
    double denom;
    arma::vec par;

    std::tie(par, denom) = nelder_mead(NM_PARAMS, inner_param.init_local, avar_wrapper,
                                       inner_param, design_info_local, distribution);

    inner_param.opt_local = to_sigmoid(par, inner_param.local_lower, inner_param.local_upper);

    return denom;
}


double asr(const arma::vec &init_coef, inner_optimization &inner_param,
           const arma::vec &glob_alloc,
           const design_info &design_info_glob, const design_info &design_info_local, int distribution) {

    // alloc: within 0 1 scale
    // local_alloc: no bound, infinity
    double nom = 0, denom = 0;
    arma::vec par, coef_sigmoid(init_coef.n_elem), local_alloc(inner_param.init_local.n_elem);

    // Nominator: global design
    if (design_info_glob.reparam) {
        coef_sigmoid = to_sigmoid(init_coef, inner_param.coef_lower, inner_param.coef_upper);
    } else {
        coef_sigmoid = init_coef;
    }

    nom = opt_crit(glob_alloc, coef_sigmoid, distribution, design_info_glob);

    if (nom == constants::BIG) {
        return - 1/constants::BIG;
    }

    // Denominator: optimal local design with the same coefficients
    inner_param.init_coef = coef_sigmoid;
    denom = local_optimal(inner_param, design_info_local, distribution);

    if (denom == constants::BIG) {
        return - 1/constants::BIG;
    }

    // Check infinity
    if (!std::isfinite(log(nom / denom))) {
        return -1/constants::BIG;
    }
    return - log(nom / denom);

}


double max_coef_asr(inner_optimization &inner_param,
                    const arma::vec &glob_alloc,
                    const design_info &design_info_glob, const design_info &design_info_local) {
    double asr_weibull, asr_lognorm, max_asr;
    arma::vec par_weibull, par_lognorm, max_par, local_opt_tmp;
    // double res_asr, max_asr = 0;
    // arma::vec par, max_par, local_opt_tmp;
    int distribution_lognorm = 2, distribution_weibull = 1;

    if (inner_param.model == 3) {

        std::tie(par_lognorm, asr_lognorm) = nelder_mead(NM_PARAMS, inner_param.init_coef, asr, inner_param,
                                                         glob_alloc, design_info_glob, design_info_local, distribution_lognorm);
        local_opt_tmp = inner_param.opt_local;

        std::tie(par_weibull, asr_weibull) = nelder_mead(NM_PARAMS, inner_param.init_coef, asr, inner_param,
                                                         glob_alloc, design_info_glob, design_info_local, distribution_weibull);



        if (asr_lognorm < asr_weibull) {
            max_asr = asr_lognorm;
            max_par = par_lognorm;
            inner_param.opt_local = local_opt_tmp;
            inner_param.opt_distribution = distribution_lognorm;
        } else {
            max_asr = asr_weibull;
            max_par = par_weibull;
            inner_param.opt_distribution = distribution_weibull;
        }


    } else if ((inner_param.model == 1) | (inner_param.model == 2)){

        std::tie(max_par, max_asr) = nelder_mead(NM_PARAMS, inner_param.init_coef, asr, inner_param,
                                             glob_alloc, design_info_glob, design_info_local, inner_param.model);

    } else {
        Rcpp::Rcout << "Model must be 1 (weibull), 2 (log normal), or 3 (model robust)." << std::endl;
        return constants::BIG;
    }

    if (design_info_glob.reparam) {
        inner_param.opt_coef = to_sigmoid(max_par, inner_param.coef_lower, inner_param.coef_upper);
    } else {
        inner_param.opt_coef = max_par;
    }

    // Evaluate and transform to asr
    if (!std::isfinite(max_asr)) {
        return -1/constants::BIG;
    }

    return max_asr;
}
