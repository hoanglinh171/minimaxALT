#ifndef COMMON_H
#define COMMON_H

#include <armadillo>
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo, RcppGSL)]]
// [[Rcpp::plugins(cpp11)]]


typedef struct {
    double maxit;
    double tol;
    double reltol;
    double alpha;
    double beta;
    double gamma;
    bool trace;
} nelder_mead_params;

namespace constants {
    const double BIG = 1e30;

    const double INTEG_WORK_ALLOC = 1000;
    const double INTEG_REL_TOL = 1e-6;
    const double INTEG_ABS_TOL = 1e-5;
    const double COND_NUM_LIM = 1e-14;
    const double UPP_P = 0.99999999999999;
    const double LOW_P = 1e-14;
    const double DERIV_TOL = 1e-7;

    const double UPP_EXP = 20;
    const double LOW_EXP = -20;
    const double LN_TOL = 1e-14;

    const double MAXIT = 500;
    const double TOL = 1e-6;
    const double RELTOL = 1e-8;
    const double ALPHA = 1;
    const double BETA = 0.5;
    const double GAMMA = 2;
    const bool TRACE = false;
}


typedef struct {
    int n_support;
    int n_factor;
    int n_unit;
    double censor_time;
    double sigma;
    double p;
    arma::vec use_cond;
    std::string opt_type;
    bool reparam;
    bool degenerate;
    double x_l;
    double x_h;
} design_info;

typedef struct {
    arma::vec opt_coef;  // sigmoid scale
    arma::vec coef_upper;
    arma::vec coef_lower;
    arma::vec init_coef;  // infinity scale
    arma::vec opt_local;  // sigmoid scale
    arma::vec local_upper;
    arma::vec local_lower;
    arma::vec init_local;  // infinity scale
    int model;
    int opt_distribution;
} inner_optimization;


/*************************************************
 * DEFINE STUCTURES OF PSO INFORMATION
 * PSO update formula:
 * w(t+1) = theta * w(t) + gamma1 * alpha1 * (pi - xi(t)) + gamma2 * alpha2 * (pg - xi(t))
 * xi(t+1) = xi(t) + w(t+1)
 * alpha1, alpha2 is uniform(0, 1)
**************************************************/

typedef struct {
    // v_step = pso_dyn.w_cur * v_step + pso_opts.c1 * arma::randu(d_swarm, n_swarm) % (p_best - swarm) +
    //          pso_opts.c2 * arma::randu(d_swarm, n_swarm) % (gb_mat - swarm);
    arma::vec var_lower;
    arma::vec var_upper;
    arma::mat init_swarm;
    int n_swarm; // 64
    int d_swarm; // 2

    int max_iter; // 100
    int early_stopping; // 10
    double tol;

    // Basic PSO Parameters
    double c1; // gamma1: 2.05
    double c2; // gamma2: 2.05
    double w0; // initial theta: 0.95
    double w1; // final theta: 0.2
    double w_var; // % iterations to update theta: 0.8
    double vk; // (max_range - min_range) / max_velocity: 2
} pso_options;




// DEFINE STUCTURES OF PSO PARAMETERS WHICH WILL CHANGE ITERATIVELY
typedef struct {
    int w_varyfor; // (int)(w_var*maxIter);
    double w_cur; // w0;
    double w_dec; // (w0 - w1)/w_varyfor;
} pso_dyn;

// DEFINE PSO RESULTS
typedef struct {
    arma::vec g_best;
    arma::vec coef_best;
    int distribution_best;
    double fg_best;
    arma::vec fg_best_hist;
    arma::mat p_best;
    arma::vec fp_best;
    arma::cube g_hist;
    arma::mat coef_best_hist;
    arma::vec distribution_best_hist;
    arma::mat model_set;
    arma::vec model_weight;
    arma::mat equivalence_data;
    double max_dd;
} pso_result;



arma::vec to_sigmoid(const arma::vec &outbound_vec, const arma::vec &lower_bound, const arma::vec &upper_bound);
arma::vec from_sigmoid(const arma::vec &inbound_vec, const arma::vec &lower_bound, const arma::vec &upper_bound);
arma::mat unique_rows(const arma::mat& m);
int compute_batch_size(int multi_start, int n_factor);

#endif // COMMON_H
