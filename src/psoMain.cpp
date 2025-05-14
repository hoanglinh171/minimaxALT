#include <psoMain.h>



double pso_obj_func(const arma::vec &particle, int design_type, inner_optimization &inner_param,
                    const design_info &design_info_glob, const design_info &design_info_local,
                    int n_multi_start, const arma::mat &init_coef_mat) {
    double obj_val = constants::BIG, feval;
    arma::vec particle_local = inner_param.opt_local, coef_col = inner_param.opt_coef;
    int opt_distribution = inner_param.opt_distribution;

    // local
    if (design_type == 1) {
        if ((inner_param.opt_distribution == 1) | (inner_param.opt_distribution == 2)) {
            obj_val = opt_crit(particle, inner_param.init_coef, inner_param.opt_distribution, design_info_glob);
        } else {
            Rcpp::Rcout << "Distribution for local optimal design must be 1 (weibull) or 2 (log normal)." << std::endl;
            return constants::BIG;
        }

        return obj_val;

    // multi-start max coef
    } else if (design_type == 2) {

        for (int j = 0; j < n_multi_start; j++) {
            inner_param.init_coef = init_coef_mat.col(j);
            feval = max_coef_asr(inner_param, particle, design_info_glob, design_info_local);

            if (feval < obj_val) {
                obj_val = feval;
                particle_local = inner_param.opt_local;
                coef_col = inner_param.opt_coef;
                opt_distribution = inner_param.opt_distribution;
            }
        }

        inner_param.opt_local = particle_local;
        inner_param.opt_coef = coef_col;
        inner_param.opt_distribution = opt_distribution;

        if (obj_val == -1/constants::BIG) {
            return constants::BIG;
        }

        return exp(-obj_val);

    } else {
        Rcpp::Rcout << "Design type must be 1 (local optimal design) or 2 (minimax design)." << std::endl;
        return constants::BIG;
    }
}


void pso_func_eval(int design_type, const arma::mat &swarm,
                   inner_optimization &original_inner_param,
                   arma::vec &f_swarm, arma::mat &coef, arma::vec &distribution_vec,
                   const design_info &design_info_glob, const design_info &design_info_local,
                   int n_multi_start, const arma::mat &init_coef_mat
                   ) {
    int n_swarm = swarm.n_cols;
    arma::vec particle_global;
    arma::mat swarm_local(3, n_swarm);
    double min_feval;

    arma::vec coef_low = original_inner_param.coef_lower;
    arma::vec coef_upp = original_inner_param.coef_upper;
    arma::vec local_low = original_inner_param.local_lower;
    arma::vec local_upp = original_inner_param.local_upper;
    arma::vec init_coef = original_inner_param.init_coef;
    arma::vec init_local = original_inner_param.init_local;
    int distribution = original_inner_param.opt_distribution;
    int model = original_inner_param.model;
    
    int batch_size = compute_batch_size(n_multi_start, design_info_glob.n_factor);

    for (int start = 0; start < n_swarm; start += batch_size) {
        Rcpp::checkUserInterrupt();  
        int end = std::min(start + batch_size, n_swarm);
    
        #pragma omp parallel for private(particle_global, min_feval) schedule(dynamic)  
        for (int i = start; i < end; i++) {
            inner_optimization inner_param =
            {init_coef, coef_upp, coef_low, init_coef, init_local, local_upp, local_low, init_local, model, distribution};

            particle_global = arma::conv_to<arma::vec>::from(swarm.col(i));
            min_feval = pso_obj_func(particle_global, design_type, inner_param, design_info_glob, design_info_local, n_multi_start, init_coef_mat);

            f_swarm(i) = min_feval;
            coef.col(i) = inner_param.opt_coef;
            swarm_local.col(i) = inner_param.opt_local;
            distribution_vec(i) = inner_param.opt_distribution;
        }
    }

}


void pso_update_dyn_para(pso_options &pso_opts, pso_dyn &pso_dyn, int iter) {
    if (iter < 0) { // INITIALIZE
        int w_varyfor = (int)(pso_opts.w_var * pso_opts.max_iter);    // number of first iterations for updating inertia weight
        pso_dyn.w_varyfor = w_varyfor;
        pso_dyn.w_cur = pso_opts.w0;
        pso_dyn.w_dec = (pso_opts.w0 - pso_opts.w1) / w_varyfor;      // Inertia weight change per iteration step
    } else { // UPDATE
        if (iter <= pso_dyn.w_varyfor) {
            pso_dyn.w_cur = pso_dyn.w_cur - pso_dyn.w_dec;
        }
    }
}


void pso_update_particle(pso_options &pso_opts, const pso_dyn &pso_dyn,
                         const arma::mat &p_best, const arma::vec &g_best,
                         arma::mat &v_step, arma::mat &swarm)
{
    int d_swarm = (int)swarm.n_rows;
    int n_swarm = (int)swarm.n_cols;

    arma::vec vel_max(d_swarm);
    vel_max = (pso_opts.var_upper - pso_opts.var_lower) / pso_opts.vk;  // max velocity allowed

    arma::mat vel_max_mat = repmat(vel_max, 1, n_swarm);
    arma::mat gb_mat = repmat(g_best, 1, n_swarm);

    GetRNGstate();

    v_step = pso_dyn.w_cur * v_step + pso_opts.c1 * arma::randu(d_swarm, n_swarm) % (p_best - swarm) +
            pso_opts.c2 * arma::randu(d_swarm, n_swarm) % (gb_mat - swarm);
    // control velocity to be smaller than max allowed
    v_step = min(v_step, vel_max_mat);
    v_step = max(v_step, (-1) * vel_max_mat);
    swarm += v_step;

    PutRNGstate();
}



void pso_check_particle(const arma::vec &var_upper, const arma::vec &var_lower, arma::mat &swarm)
{
    int d_swarm = (int)swarm.n_rows;
    int n_swarm = (int)swarm.n_cols;

    arma::mat var_upp_mat = repmat(var_upper, 1, n_swarm);
    arma::mat var_low_mat = repmat(var_lower, 1, n_swarm);

    arma::mat swarm_tmp;
    // UPDATE POSITION

    GetRNGstate();

    arma::umat swarm_change;
    swarm_change = find(swarm > var_upp_mat);
    swarm_tmp = arma::pow(arma::randu(d_swarm, n_swarm), 5e-2) % (var_upp_mat - var_low_mat) + var_low_mat;
    for (arma::uword i = 0; i < swarm_tmp.n_cols; i++) {
        for (arma::uword j = 0; j < swarm_tmp.n_rows; j++) {
            double coin = arma::randu<double>();
            if (coin > 0.5) {
                swarm_tmp(j, i) = var_upp_mat(j, i);
            }
        }
    }
    swarm.elem(swarm_change) = swarm_tmp.elem(swarm_change);

    swarm_change = find(swarm < var_low_mat);
    swarm_tmp = (1.0 - arma::pow(1.0 - arma::randu(d_swarm, n_swarm), 5e-2)) % (var_upp_mat - var_low_mat) + var_low_mat;
    for (arma::uword i = 0; i < swarm_tmp.n_cols; i++) {
        for (arma::uword j = 0; j < swarm_tmp.n_rows; j++) {
            double coin = arma::randu<double>();
            if (coin > 0.5) {
                swarm_tmp(j, i) = var_low_mat(j, i);
            }
        }
    }
    swarm.elem(swarm_change) = swarm_tmp.elem(swarm_change);

    PutRNGstate();

}


void pso_main(int design_type, pso_options &pso_opts, inner_optimization &inner_param,
    design_info& design_info_local, design_info& design_info_glob, pso_result &pso_result,
    int n_multi_start, const arma::mat &init_coef_mat, bool verbose)
{
    // GET PSO PARAMETERS
    int n_swarm = pso_opts.n_swarm;
    int d_swarm = pso_opts.d_swarm;  // number of elements of optimizing variable vector
    int max_iter = pso_opts.max_iter;
    int early_stopping = pso_opts.early_stopping;
    double tol = pso_opts.tol;
    arma::vec var_lower  = pso_opts.var_lower;
    arma::vec var_upper  = pso_opts.var_upper;

    // DECLARE VARIABLES
    arma::mat swarm(d_swarm, n_swarm), v_step(d_swarm, n_swarm), p_best(d_swarm, n_swarm);//, GrBest(nGroup, d_swarm);
    arma::vec vel_max(d_swarm);
    arma::mat coef(inner_param.coef_upper.n_elem, n_swarm);
    arma::vec g_best(d_swarm), coef_best(inner_param.coef_upper.n_elem);
    arma::vec f_swarm(n_swarm), distribution_swarm(n_swarm), fp_best(n_swarm);//, fGrBest(nGroup);
    double fg_best;
    int distribution_best;
    arma::uword g_best_idx;
    arma::vec fg_best_hist(max_iter + 1, arma::fill::zeros);
    arma::mat coef_best_hist(inner_param.coef_upper.n_elem, max_iter + 1);
    arma::vec distribution_best_hist(max_iter + 1, arma::fill::zeros);
    arma::cube g_hist(d_swarm, n_swarm, max_iter + 1);
    pso_dyn pso_dyn;


    Rcpp::Rcout << "PSO Loop: Initializing..    " << std::endl;

    // INITIALIZE RANDOM SWARM
    swarm = pso_opts.init_swarm;

    // INITIALIZE VELOCITY
    v_step.fill(0);

    // INITIALIZE OBJECTIVE FUNCTION VALUES (CHANGE f_swarm, local_swarm, COEF)
    pso_func_eval(design_type, swarm, inner_param, f_swarm, coef, distribution_swarm, design_info_glob, design_info_local,
                  n_multi_start, init_coef_mat);

    // INITIALIZE LOCAL BEST
    fp_best = f_swarm;	p_best = swarm;

    // INITIALIZE GLOBAL BEST
    fg_best = fp_best.min(g_best_idx);
    g_best = p_best.col(g_best_idx);
    coef_best = coef.col(g_best_idx);
    distribution_best = distribution_swarm(g_best_idx);

    // INITIALIZE PSO DYNAMIC PARAMETERS
    pso_update_dyn_para(pso_opts, pso_dyn, -1);

    // SAVE INITIAL GLOBAL BEST VALUE
    fg_best_hist(0) = fg_best;
    coef_best_hist.col(0) = coef_best;
    g_hist.slice(0) = swarm;
    distribution_best_hist(0) = distribution_best;

    Rcpp::Rcout << "OK" << std::endl;
    Rcpp::Rcout << std::endl;
    

    /* -- START PSO LOOP -- */
    int t = 0;
    int k = 0;

    try {

        while (t < max_iter) {
            // Rcpp::checkUserInterrupt();

            auto start = std::chrono::high_resolution_clock::now();

            if (t == 0)  Rcpp::Rcout << "PSO Loop: Updating ..    " << std::endl;
            if (verbose) {
                Rcpp::Rcout << "Processing: " << t + 1 << "/" << max_iter << std::endl;
                Rcpp::Rcout << "-------------------------------------------------------------------------------" << std::endl;
            }
            if (t == (max_iter - 1)) Rcpp::Rcout << std::endl;

            // UPDATE VELOCITY
            pso_update_particle(pso_opts, pso_dyn, p_best, g_best, v_step, swarm);
            pso_check_particle(var_upper, var_lower, swarm);

            // UPDATE OBJECTIVE FUNCTION VALUES
            pso_func_eval(design_type, swarm, inner_param, f_swarm, coef, distribution_swarm, design_info_glob, design_info_local,
                        n_multi_start, init_coef_mat);

            // UPDATE THE LOCAL AND GLOBAL BEST
            if (any(f_swarm < fp_best)) {
                arma::uvec colchange = find(f_swarm < fp_best);
                fp_best.elem(colchange) = f_swarm.elem(colchange);
                p_best.cols(colchange) = swarm.cols(colchange);
            }
            if (min(fp_best) < fg_best) {
                fg_best = fp_best.min(g_best_idx);
                g_best = p_best.col(g_best_idx);
                coef_best = coef.col(g_best_idx);  //pso_dyn.succ_GB = 1;
                distribution_best = distribution_swarm(g_best_idx);
            }

            // UPDATE PSO DYNAMIC PARAMETERS
            pso_update_dyn_para(pso_opts, pso_dyn, t);

            // SAVE CURRENT GLOBAL BEST VALUE
            fg_best_hist(t+1) = fg_best;
            coef_best_hist.col(t+1) = coef_best;
            distribution_best_hist(t+1) = distribution_best;
            g_hist.slice(t+1) = swarm;


            // Only show progress if design_type = 2


            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            if (verbose) Rcpp::Rcout << "Time taken by the iteration: " 
                                    << duration.count() / 1000 << " seconds" << std::endl;
            if (verbose) Rcpp::Rcout << std::endl;

                // CHECK STOPPING CRITERION

                if ((t + 1 > 0) & ((t + 1) % early_stopping == 0)) {

                    // Assign to struct
                    pso_result.g_best = g_best;
                    pso_result.coef_best = coef_best;
                    pso_result.distribution_best = distribution_best;
                    pso_result.fg_best = fg_best;
                    pso_result.fg_best_hist = fg_best_hist.subvec(0, t + 1);
                    pso_result.p_best = p_best;
                    pso_result.fp_best = fp_best;
                    pso_result.g_hist = g_hist.slices(0, t + 1);
                    pso_result.coef_best_hist = coef_best_hist.cols(0, t + 1);
                    pso_result.distribution_best_hist = distribution_best_hist.subvec(0, t + 1);

                    // Check equivalence theorem
                    arma::mat model_set (t + 2, coef_best.n_elem + 1);
                    model_set.submat(0, 0, model_set.n_rows - 1, model_set.n_cols - 2) = coef_best_hist.cols(0, t + 1).t();
                    model_set.col(coef_best.n_elem) = distribution_best_hist.subvec(0, t + 1);

                    model_set = unique_rows(model_set);
                    pso_result.model_set = model_set;
                    equivalence_plot_data(g_best, design_info_glob, pso_result);

                    if(verbose) Rcpp::Rcout << "###### Max directional derivative: " << pso_result.max_dd << " ######" << std::endl;

                    // If the design is optimal, stop
                    if (abs(pso_result.max_dd - 1) < tol) {
                        t = max_iter;
                    }
                }


                if (verbose) {
                    Rcpp::Rcout << "Objective value: " << fg_best << std::endl;
                    Rcpp::Rcout << "Levels: " << g_best.subvec(0, design_info_glob.n_factor * design_info_glob.n_support - 1).t() << std::endl;
                    Rcpp::Rcout << "Proportion: " << softmax(g_best.subvec(design_info_glob.n_factor * design_info_glob.n_support, g_best.n_elem - 1)).t() << std::endl;
                    Rcpp::Rcout << "Distribution: " << distribution_best << std::endl;
                    Rcpp::Rcout << "Parameters: " << coef_best.t() << std::endl;

                    Rcpp::Rcout << std::endl;
                }

                t++;
                k++;

            }

        /* -- FINISH PSO LOOP -- */
    } catch (...) {
        Rcpp::Rcout << "###### Interrupted by user ######" << std::endl;
        Rcpp::Rcout << "Objective value: " << fg_best << std::endl;
        Rcpp::Rcout << "Levels: " << g_best.subvec(0, design_info_glob.n_factor * design_info_glob.n_support - 1).t() << std::endl;
        Rcpp::Rcout << "Proportion: " << softmax(g_best.subvec(design_info_glob.n_factor * design_info_glob.n_support, g_best.n_elem - 1)).t() << std::endl;
        Rcpp::Rcout << "Distribution: " << distribution_best << std::endl;
        Rcpp::Rcout << "Parameters: " << coef_best.t() << std::endl;
    }

    /* -- OUTPUT -- */
    pso_result.g_best = g_best;
    pso_result.coef_best = coef_best;
    pso_result.distribution_best = distribution_best;
    pso_result.fg_best = fg_best;
    pso_result.fg_best_hist = fg_best_hist.subvec(0, k);
    pso_result.p_best = p_best;
    pso_result.fp_best = fp_best;
    pso_result.g_hist = g_hist.slices(0, k);
    pso_result.coef_best_hist = coef_best_hist.cols(0, k);
    pso_result.distribution_best_hist = distribution_best_hist.subvec(0, k);

    // Check equivalence theorem
    arma::mat model_set (k + 1, coef_best.n_elem + 1);
    model_set.submat(0, 0, model_set.n_rows - 1, model_set.n_cols - 2) = coef_best_hist.cols(0, k).t();
    model_set.col(coef_best.n_elem) = distribution_best_hist.subvec(0, k);

    model_set = unique_rows(model_set);
    pso_result.model_set = model_set;
    equivalence_plot_data(g_best, design_info_glob, pso_result);

    /* -- END -- */
}
