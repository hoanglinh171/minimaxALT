test_that("locally optimal design is reproducible with the same seed", {
    
    result1 <- find_optimal_alt(
        design_type = "locally",
        distribution = "weibull",
        design_info = design_info_test,
        pso_info = pso_info_test,
        coef = coef_reparam,
        verbose = FALSE,
        seed = 123
    )
    
    result2 <- find_optimal_alt(
        design_type = "locally",
        distribution = "weibull",
        design_info = design_info_test,
        pso_info = pso_info_test,
        coef = coef_reparam,
        verbose = FALSE,
        seed = 123
    )
    
    expect_identical(result1, result2)
})


test_that("locally optimal design with two factors is reproducible with the same seed", {
    
    result1 <- find_optimal_alt(
        design_type = "locally",
        distribution = "weibull",
        design_info = design_info_2f_test,
        pso_info = pso_info_test,
        coef = coef_2f_reparam,
        verbose = FALSE,
        seed = 123
    )
    
    result2 <- find_optimal_alt(
        design_type = "locally",
        distribution = "weibull",
        design_info = design_info_2f_test,
        pso_info = pso_info_test,
        coef = coef_2f_reparam,
        verbose = FALSE,
        seed = 123
    )
    
    expect_identical(result1, result2)
})



test_that("minimax design is reproducible with the same seed", {
    
    result1 <- find_optimal_alt(
        design_type = "minimax",
        distribution = "lognormal",
        design_info = design_info_test,
        pso_info = pso_info_test,
        coef_upper = coef_upper,
        coef_lower = coef_lower,
        verbose = FALSE,
        seed = 123
    )
    
    result2 <- find_optimal_alt(
        design_type = "minimax",
        distribution = "lognormal",
        design_info = design_info_test,
        pso_info = pso_info_test,
        coef_upper = coef_upper,
        coef_lower = coef_lower,
        verbose = FALSE,
        seed = 123
    )
    
    expect_identical(result1, result2)
})


test_that("optimality checkk is reproducible with the same seed", {
    
    result1 <- opt_check <- check_optimality(
        best_design = best_design_test,
        model_set = model_set_test,
        design_info = design_info_test,
        seed = 123
    )
    
    result2 <- opt_check <- check_optimality(
        best_design = best_design_test,
        model_set = model_set_test,
        design_info = design_info_test,
        seed = 123
    )
    
    expect_identical(result1, result2)
})

