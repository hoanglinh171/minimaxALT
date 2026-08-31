test_that("locally optimal design can be applied for two-factor design", {
    
    design_result <- find_optimal_alt(
        design_type = "locally", 
        distribution = "weibull", 
        design_info = design_info_2f_test, 
        pso_info = pso_info_test,
        coef = coef_2f_reparam,
        verbose = FALSE
    )
    
    extracted_2f_design <- extract_design(design_result)
    
    expect_s3_class(design_result, "OptimalALT")
    expect_equal(ncol(extracted_2f_design), design_info_2f_test$n_support)
    expect_equal(nrow(extracted_2f_design), design_info_2f_test$n_factor + 1)
    
})


test_that("minimax design can be applied for two-factor design", {
    
    design_result <- find_optimal_alt(
        design_type = "minimax", 
        distribution = "lognormal", 
        design_info = design_info_2f_test, 
        pso_info = pso_info_test,
        coef_lower = c(10^-6, 10^-6, 0.7),
        coef_upper = c(10^-4, 10^-4, 0.9),
        highest_level = TRUE,
        verbose = FALSE
    )
    
    extracted_2f_design <- extract_design(design_result)
    
    expect_s3_class(design_result, "OptimalALT")
    expect_equal(ncol(extracted_2f_design), design_info_2f_test$n_support)
    expect_equal(nrow(extracted_2f_design), design_info_2f_test$n_factor + 1)
    
})