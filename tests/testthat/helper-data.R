n_swarm = 3
d_swarm = 5
x_l = 0
x_h = 1
coef_lower = c(10^-6, 0.7)
coef_upper = c(10^-3, 0.9)
coef_reparam = c(0.001, 0.9)
coef = c(9, -0.1)
coef_2f = c(9, -0.1, -0.05)
coef_2f_reparam = c(0.001, 0.1, 1)

design_info_test <- set_design_info(
    k_levels=3, 
    j_factor=1, 
    n_unit=300,
    censor_time=183, 
    p=0.1, 
    use_cond=0, 
    sigma=0.6
)

design_info_original_test <- design_info_test
design_info_original_test$reparam <- FALSE

design_info_2f_test <- set_design_info(
    k_levels=3, 
    j_factor=2, 
    n_unit=300,
    censor_time=183, 
    p=0.1, 
    use_cond=c(0, 0), 
    sigma=0.6
)

pso_info_test <- pso_setting(
    n_swarm=n_swarm, 
    max_iter=10, 
    early_stopping=5, 
    tol=0.0001
)

init_values_test <- initialize_values(
    init_swarm = matrix(runif(n_swarm * d_swarm, min = x_l, max = x_h), 
                        nrow = n_swarm),
    init_local = c(1, 0.7, 0.3),
    init_coef_mat = t(replicate(4, runif(2, min = coef_lower, max = coef_upper)))
)

design_result_test <- find_optimal_alt(
    design_type = "minimax", 
    distribution = "lognormal", 
    design_info = design_info_test, 
    pso_info = pso_info_test, 
    coef_lower = coef_lower,
    coef_upper = coef_upper,
    highest_level = TRUE,
    verbose = FALSE,
    n_threads = 4
)


best_design_test <- rbind(
    c(1, 0.7, 0.3),
    c(0.2, 0.1, 0.7)
)

model_set_test <- rbind(
    c(10^-6, 0.7, "weibull"),
    c(10^-3, 0.9, "weibull")
)


