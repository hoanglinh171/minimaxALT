devtools::load_all()
library(minimaxALT)


init_coef_mat = rbind(
  c(1.0005e-4, 0.7000001),
  c(1.0005e-4, 0.8999999),
  c(0.999999e-1, 0.8999999),
  c(0.999999e-1, 0.7000001)
)


init_values = initialize_values(init_coef_mat = init_coef_mat)


design_info <- set_design_info(k_levels=4, j_factor=1, n_unit=300,
                               censor_time=183, p=0.1, use_cond=c(0),
                               sigma=0.6, x_l = 0.1, x_h = 1)

pso_info <- pso_setting(n_swarm=64, max_iter=180, early_stopping=10, tol=0.01)


res <- find_optimal_alt(design_type=2, distribution=3, design_info=design_info,
                                   pso_info=pso_info,
                                   # coef=c(1.82*10^-6, 3.17*10^-6, 1),
                                   coef_lower = c(1e-4, 0.7),
                                   coef_upper = c(1e-1, 0.9),
                                   # highest_level = FALSE,
                                   init_values = init_values,
                                   verbose = TRUE, n_threads = 16
)


plot(res, x_l = 0.1, x_h =1)

summary(res)

### ---------------------------------------------------------------------------
### Equivalence theorem


equi_check <- check_equivalence_theorem (best_particle = as.vector(res$g_best), model_set = res$model_set,
                                         design_info = design_info)


equi_check$max_directional_derivative
equi_check$model_weight
