

# init_coef <- rbind(
#   c(1e-2, 1e-2, 1),
#   c(1e-6, 5e-3, 1),
#   c(1e-2, 1e-2, 0.7),
#   c(1e-6, 5e-3, 0.7),
#   c(1e-6, 1e-6, 1),
#   c(1e-6, 1e-6, 0.7)
# )
# 
# init_values = initialize_values(init_coef_mat = init_coef)

design_info <- set_design_info(k_levels=3, j_factor=1, n_unit=300,
                               censor_time=170, p=0.1, use_cond=c(0),
                               sigma=0.6, x_l = 0, x_h = 1)

pso_info <- pso_setting(n_swarm=64, max_iter=128, early_stopping=10, tol=0.01)


res <- find_optimal_alt(design_type=1, distribution=1, design_info=design_info,
                                   pso_info=pso_info,
                                   coef=c(1.82*10^-6, 1),
                                   # coef_lower = c(1e-6, 0.7),
                                   # coef_upper = c(1e-2, 1),
                                   # highest_level = FALSE,
                                   # init_values = init_values,
                                   verbose = TRUE, n_threads = 1
)


plot(res, x_l = 0, x_h =1, nlevels = 10)

summary(res)

### ---------------------------------------------------------------------------
### Equivalence theorem


equi_check <- check_equivalence_theorem (best_particle = as.vector(res$g_best), model_set = res$model_set,
                                         design_info = design_info)


equi_check$max_directional_derivative
equi_check$model_weight
