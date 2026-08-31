test_that("design info returns a DesignInfo object", {
    
    design_info <- set_design_info(
        k_levels=3, 
        j_factor=1, 
        n_unit=300,
        censor_time=183, 
        p=0.1, 
        use_cond=0, 
        sigma=0.6
    )
    
    expect_s3_class(design_info, "DesignInfo")
    
})

test_that("pso info returns a PSOInfo object", {
    
    pso_info <- pso_setting(
        n_swarm=n_swarm, 
        max_iter=30, 
        early_stopping=10, 
        tol=0.0001
    )
    
    expect_s3_class(pso_info, "PSOInfo")

})


test_that("initialize values returns a InitialValue object", {
    
    init_values <- initialize_values(
        init_swarm = matrix(runif(n_swarm * d_swarm, min = x_l, max = x_h), 
                            nrow = n_swarm),
        init_local = c(1, 0.7, 0.3),
        init_coef_mat = t(replicate(4, runif(2, min = coef_lower, max = coef_upper)))
    )
    
    expect_s3_class(init_values, "InitialValue")
    
})


test_that("optimal ALT design returns an OptimalALT object", {
    
    design_result <- find_optimal_alt(
        design_type = "locally", 
        distribution = "lognormal", 
        design_info = design_info_test, 
        pso_info = pso_info_test,
        coef = coef_reparam,
        highest_level = TRUE,
        verbose = FALSE
    )
    
    expect_s3_class(design_result, "OptimalALT")

})


test_that("check optimality returns an OptimalCheck object", {
    
    opt_check <- check_optimality(
        best_design = best_design_test,
        model_set = model_set_test,
        design_info = design_info_test
    )
    
    expect_s3_class(opt_check, "OptimalityCheck")
    
})

