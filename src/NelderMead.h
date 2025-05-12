#ifndef NELDERMEAD_H
#define NELDERMEAD_H

#include "common.h"
#include "asymptoticVariance.h"

const nelder_mead_params NM_PARAMS = {
    constants::MAXIT,
    constants::TOL,
    constants::RELTOL,
    constants::ALPHA,
    constants::BETA,
    constants::GAMMA,
    constants::TRACE
};

template <typename Func, typename Param, typename... Args>
std::tuple<arma::vec, double> nelder_mead(const nelder_mead_params &nm_params,
                                          const Param &init, const Func &func, Args&... args) {

    double maxit = nm_params.maxit;
    double tol = nm_params.tol;
    double reltol = nm_params.reltol;
    double alpha = nm_params.alpha;
    double beta = nm_params.beta;
    double gamma = nm_params.gamma;
    bool trace = nm_params.trace;

    int n = init.n_elem;
    double step, trystep;
    arma::mat simplex(n, n+1);
    arma::vec fvec(n+1);
    arma::vec centroid(n), r(n), e(n), c(n);
    double fr = 0, fe = 0, fc = 0;
    arma::uvec order;


    // initialize simplex
    step = 0.1 * abs(init).max();
    if (trace) std::cout << "step: " << step << std::endl;

    if (step == 0) step = 0.1;
    simplex.col(0) = init;  // let the first row to be input initial
    for (int i = 1; i <= n; ++i) {
        simplex.col(i) = init;

        trystep = step;
        while (simplex(i-1, i) == init(i-1)) {
            simplex(i-1, i) = init(i-1) + trystep;
            trystep *= 10;
        }
    }

    // simplex.print("Build simplex: ");

    // calculate obj value and sort
    for (arma::uword j = 0; j < simplex.n_cols; j++) {
        fvec(j) = func(simplex.col(j), args...);
    }

    double convrel = reltol * (fabs(fvec(0)) + reltol);
    // double convrel = reltol;

    if (trace) std::cout << "convrel: " << convrel << std::endl;
    if (trace) fvec.t().print("Build: ");



    // ------------------- START LOOP -----------------------------------------------------
    int iter = 0;
    while (iter < maxit) {
        if (fvec.has_nan()) {
            simplex.raw_print("Simplex: ");
            fvec.t().raw_print("Fvec: ");
        }
        order = arma::sort_index(fvec);
        simplex = simplex.cols(order);
        fvec = fvec(order);

        if (fvec(n) <= fvec(0) + convrel || fabs(fvec(0)) <= tol) break;

        centroid = arma::mean(simplex.cols(0, n-1), 1);

        // Reflection
        r = (1 + alpha) * centroid - alpha * simplex.col(n);
        arma::vec R = arma::conv_to<arma::vec>::from(r);
        fr = func(R, args...);

        if (fr < fvec(0)) {

            // Extend
            e = gamma * r + (1 - gamma) * centroid;
            arma::vec E = arma::conv_to<arma::vec>::from(e);
            fe = func(E, args...);

            if (fe < fr) {
                // Replace w with e
                simplex.col(n) = e;
                fvec(n) = fe;

                if (trace) fvec.t().print("Extension: ");

            } else {
                // Replace w with r
                simplex.col(n) = r;
                fvec(n) = fr;

                if (trace) fvec.t().print("Reflection: ");

            }
        } else if (fr < fvec(n)) {
            // Replace w with r
            simplex.col(n) = r;
            fvec(n) = fr;

            if (trace) fvec.t().print("Hi Reduction");

        } else {
            // Contraction
            c = (1 - beta) * simplex.col(n) + beta * centroid;
            arma::vec C = arma::conv_to<arma::vec>::from(c);
            fc = func(C, args...);

            if (fc < fvec(n)) {
                // Replace w with c
                simplex.col(n) = c;
                fvec(n) = fc;

                if (trace) fvec.t().print("Lo Reduction");

            } else {
                // Shrink
                for (arma::uword k = 1; k < fvec.n_elem; k++) {
                    simplex.col(k) = beta * (simplex.col(k) - simplex.col(0)) + simplex.col(0);
                    fvec(k) = func(simplex.col(k), args...);
                }

                if (trace) fvec.t().print("Shrink: ");
            }
        }
        iter ++;
    }

    if (trace) std::cout << "Iteration: " << iter << std::endl;

    return std::make_tuple(simplex.col(0), fvec(0));
}


double avar_wrapper(const arma::vec &init_local, inner_optimization &inner_param,
                    const design_info &design_info_local, int distribution);

double local_optimal(inner_optimization &inner_param, const design_info &design_info_local, int distribution);

double asr(const arma::vec &init_coef, inner_optimization &inner_param,
           const arma::vec &glob_alloc,
           const design_info &design_info_glob, const design_info &design_info_local,
           int distribution);

double max_coef_asr(inner_optimization &inner_param,
                    const arma::vec &glob_alloc,
                    const design_info &design_info_glob, const design_info &design_info_local);

#endif // NELDERMEAD_H
