#ifndef PSO_H
#define PSO_H

#include "common.h"
#include "NelderMead.h"
#include "asymptoticVariance.h"
#include "equivalenceTheorem.h"


double pso_obj_func(const arma::vec &particle, int design_type, inner_optimization &inner_param,
                    const design_info &design_info_glob, const design_info &design_info_local,
                    int n_multi_start, const arma::mat &init_coef_mat);

void pso_func_eval(int design_type, const arma::mat &swarm,
                   inner_optimization &original_inner_param,
                   arma::vec &f_swarm, arma::mat &coef, arma::vec &distribution_vec,
                   const design_info &design_info_glob, const design_info &design_info_local,
                   int n_multi_start, const arma::mat &init_coef_mat
                   );

void pso_update_dyn_para(pso_options &pso_opts, pso_dyn &pso_dyn, int iter);

void pso_update_particle(pso_options &pso_opts, const pso_dyn &pso_dyn,
                         const arma::mat &p_best, const arma::vec &g_best,
                         arma::mat &v_step, arma::mat &swarm);

void pso_check_particle(const arma::vec &var_upper, const arma::vec &var_lower, arma::mat &swarm);

void pso_main(int design_type, pso_options &pso_opts, inner_optimization &inner_param,
              design_info& design_info_local, design_info& design_info_glob, pso_result &pso_result,
              int n_multi_start, const arma::mat &init_coef_mat, bool verbose);


#endif // PSO_H
