design_info <- set_design_info(k_levels=3, j_factor=1, n_unit=300, censor_time=183, 
                               p=0.1, use_cond=c(0), sigma=0.6, x_l = 0, x_h = 1)

pso_info <- pso_setting(n_swarm=20, max_iter=100, early_stopping=10, tol=0.0001)

res <- find_optimal_alt(design_type=2, distribution=1, 
                        design_info=design_info, 
                        pso_info=pso_info, 
                        coef = c(1e-3, 0.9),
                        verbose = TRUE, n_threads = 10
                        )

summary(res)
plot(res)

