design_info <- set_design_info(k_levels=3, j_factor=2, n_unit=170, 
                               censor_time=1000, p=0.001, use_cond=c(0, 0), sigma=0.6734, 
                               reparam=FALSE)

pso_info <- pso_setting(n_swarm=16, max_iter=30, early_stopping=10, 
                        tol=0.0001)

set.seed(5)
example2_locally_optimal <- find_optimal_alt(
    design_type="locally", distribution="weibull", 
    design_info=design_info, pso_info=pso_info, 
    coef=c(15.808, -11.249, -3.006), verbose = FALSE)


summary(example2_locally_optimal)

plot(example2_locally_optimal, x_l=0, x_h=1)