design_info <- set_design_info(k_levels=3, j_factor=1, n_unit=300,
                               censor_time=183, p=0.1, use_cond=c(0),
                               sigma=0.6, x_l = 0.1, x_h = 1)

pso_info <- pso_setting(n_swarm=90, max_iter=180, early_stopping=5, tol=0.0001)

res <- find_optimal_alt(design_type=1, distribution=1, design_info=design_info,
                        pso_info=pso_info,
                        coef=c(1.82*10^-6, 1),
                        # coef_lower = c(1e-5, 0.7),
                        # coef_upper = c(1e-4, 0.9),
                        # highest_level = FALSE,
                        # init_values = init_values,
                        verbose = TRUE, n_threads = 10
                        )

plot(res)

