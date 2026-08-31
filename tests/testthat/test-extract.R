test_that("extract_design returns a matrix of levels and weights", {
    
    extracted_design <- extract_design(design_result_test)
    
    expect_equal(ncol(extracted_design), design_info_test$n_support)
    expect_equal(nrow(extracted_design), design_info_test$n_factor + 1)
    expect_equal(sum(extracted_design["W", ]), 1, tolerance = 1e-10)
    
})

