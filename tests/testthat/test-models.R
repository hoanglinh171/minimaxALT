test_that("invalid distribution values are rejected", {
    
    expect_error(
        find_optimal_alt(
            design_type = "locally",
            distribution = "normal",
            design_info = design_info_test,
            pso_info = pso_info_test,
            coef = coef_reparam,
            verbose = FALSE
        )
    )
    
    expect_error(
        find_optimal_alt(
            design_type = "locally",
            distribution = c("weibull", "lognormal"),
            design_info = design_info_test,
            pso_info = pso_info_test,
            coef = coef_reparam,
            verbose = FALSE
        )
    )
    
    expect_error(
        find_optimal_alt(
            design_type = "locally",
            distribution = NA_character_,
            design_info = design_info_test,
            pso_info = pso_info_test,
            coef = coef_reparam,
            verbose = FALSE
        )
    )
    
    expect_error(
        find_optimal_alt(
            design_type = "locally",
            distribution = "",
            design_info = design_info_test,
            pso_info = pso_info_test,
            coef = coef_reparam,
            verbose = FALSE
        )
    )
})

test_that("original coefficients are not accepted with reparam = TRUE", {
    
    expect_error(
        minimax_alt(
            design_type = "locally",
            distribution = "weibull",
            design_info = design_info_test,
            pso_info = pso_info_test,
            coef = coef,
            verbose = FALSE
        )
    )
})


test_that("both original and reparameterized coefficients are accepted with reparam = FALSE", {
    
    expect_no_error(
        find_optimal_alt(
            design_type = "locally",
            distribution = "weibull",
            design_info = design_info_original_test,
            pso_info = pso_info_test,
            coef = coef,
            verbose = FALSE
        )
    )
    
    expect_no_error(
        find_optimal_alt(
            design_type = "locally",
            distribution = "weibull",
            design_info = design_info_original_test,
            pso_info = pso_info_test,
            coef = coef_reparam,
            verbose = FALSE
        )
    )
})


