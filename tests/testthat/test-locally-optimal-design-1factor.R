design_info <- set_design_info(k_levels=2, j_factor=1, n_unit=300, 
                               censor_time=183, p=0.1, use_cond=0, sigma=0.6)

pso_info <- pso_setting(n_swarm=16, max_iter=30, early_stopping=10, 
                        tol=0.0001)

set.seed(4)
example1_locally_optimal <- find_optimal_alt(
    design_type="locally", distribution="weibull", 
    design_info=design_info, pso_info=pso_info, 
    coef = c(0.001, 0.9),
    verbose = FALSE
)

summary(example1_locally_optimal)

plot(example1_locally_optimal, x_l=0, x_h=1)