#### Transformation functions --------------------------------------------------

## Transform coefficients by sigmoid function

get_outbound_sigmoid <- function(coef_vec, coef_lower, coef_upper) {
    
  stopifnot(is.vector(coef_vec),
            all(is.finite(coef_vec)),
            all(is.finite(coef_lower)),
            all(is.finite(coef_upper)),
            (length(coef_vec) == length(coef_lower)),
            all(coef_vec >= coef_lower),
            all(coef_vec <= coef_upper)
            )
  
  return(transform_sigmoid(coef_vec, coef_lower, coef_upper))
}

## Get proportion of designs

get_proportion <- function(transform_prop) {
    
  stopifnot(is.vector(transform_prop),
            all(is.finite(transform_prop)),
            all(transform_prop >= 0),
            all(transform_prop <= 1)
  )
  
  return(transform_proportion(transform_prop))
}

## Get transformed proportion for PSO

get_transform_prop <- function(prop) {
  stopifnot(is.vector(prop),
            all(is.finite(prop)),
            all(prop >= 0),
            all(prop <= 1),
            sum(prop) == 1
  )
  
  transform_prop = prop[1:(length(prop) - 1)]
  if (length(transform_prop) > 1) {
    denom = 1
    for (i in 2:length(transform_prop)) {
      denom = denom - prop[i-1]
      transform_prop[i] = prop[i] / denom 
    }
  }
  
  return(transform_prop)
}


#### Validate inputs -----------------------------------------------------------

check_design_info <- function(design_info) {
    
    if (!inherits(design_info, "DesignInfo")) {
        stop("`design_info` must be an object of class `DesignInfo`.")
    }
    
    # k_levels
    if (!is.numeric(design_info$n_support) || 
        length(design_info$n_support) != 1 ||
        !is.finite(design_info$n_support) || 
        design_info$n_support < 1 || 
        design_info$n_support != as.integer(design_info$n_support)) {
        stop("`k_levels` must be a positive integer.")
    }
    
    # j_factor
    if (!is.numeric(design_info$n_factor) || 
        length(design_info$n_factor) != 1 ||
        !is.finite(design_info$n_factor) || 
        design_info$n_factor < 1 || 
        design_info$n_factor != as.integer(design_info$n_factor)) {
        stop("`j_factor` must be a positive integer.")
    }
    
    # n_unit
    if (!is.numeric(design_info$n_unit) || length(design_info$n_unit) != 1 ||
        !is.finite(design_info$n_unit) || design_info$n_unit < 1 || 
        design_info$n_unit != as.integer(design_info$n_unit)) {
        stop("`n_unit` must be a positive integer.")
    }
    
    # censor_time
    if (!is.numeric(design_info$censor_time) || 
        length(design_info$censor_time) != 1 ||
        !is.finite(design_info$censor_time) || 
        design_info$censor_time <= 0) {
        stop("`censor_time` must be a positive numeric value.")
    }
    
    # p
    if (!is.numeric(design_info$p) || length(design_info$p) != 1 ||
        !is.finite(design_info$p) || design_info$p < 0 || design_info$p > 1) {
        stop("`p` must be a numeric value between 0 and 1.")
    }
    
    # use_cond
    if (!is.numeric(design_info$use_cond) || 
        any(!is.finite(design_info$use_cond))) {
        stop("`use_cond` must be a numeric vector containing only finite values.")
    }
    
    if (length(design_info$use_cond) != design_info$n_factor) {
        stop("The length of `use_cond` must equal `j_factor`.")
    }
    
    if (any(design_info$use_cond < 0 | design_info$use_cond > 1)) {
        stop("All values of `use_cond` must be between 0 and 1.")
    }
    
    # sigma
    if (!is.numeric(design_info$sigma) || length(design_info$sigma) != 1 ||
        !is.finite(design_info$sigma) || design_info$sigma <= 0) {
        stop("`sigma` must be a positive numeric value.")
    }
    
    # x_l
    if (!is.numeric(design_info$x_l) || length(design_info$x_l) != 1 ||
        !is.finite(design_info$x_l) || 
        design_info$x_l < 0 || design_info$x_l > 1) {
        stop("`x_l` must be a single finite numeric value between 0 and 1.")
    }
    
    # x_h
    if (!is.numeric(design_info$x_h) || length(design_info$x_h) != 1 ||
        !is.finite(design_info$x_h) || 
        design_info$x_h < 0 || design_info$x_h > 1) {
        stop("`x_h` must be a single finite numeric value between 0 and 1.")
    }
    
    # stress range
    if (design_info$x_l > design_info$x_h) {
        stop("`x_l` must be smaller than `x_h`.")
    }
    
    # reparam
    if (!is.logical(design_info$reparam) || length(design_info$reparam) != 1 ||
        is.na(design_info$reparam)) {
        stop("`reparam` must be a single TRUE or FALSE value.")
    }
    
    invisible(TRUE)
    
}


check_pso_info <- function(pso_info) {
    
    if (!inherits(pso_info, "PSOInfo")) {
        stop("`pso_info` must be an object of class `PSOInfo`.")
    }
    
    # Number of particles
    if (!is.numeric(pso_info$n_swarm) ||
        length(pso_info$n_swarm) != 1 ||
        !is.finite(pso_info$n_swarm) ||
        pso_info$n_swarm < 1 ||
        pso_info$n_swarm != as.integer(pso_info$n_swarm)) {
        stop("`n_swarm` must be a positive integer.")
    }
    
    # Maximum number of iterations
    if (!is.numeric(pso_info$max_iter) ||
        length(pso_info$max_iter) != 1 ||
        !is.finite(pso_info$max_iter) ||
        pso_info$max_iter < 1 ||
        pso_info$max_iter != as.integer(pso_info$max_iter)) {
        stop("`max_iter` must be a positive integer.")
    }
    
    # Early stopping frequency
    if (!is.numeric(pso_info$early_stopping) ||
        length(pso_info$early_stopping) != 1 ||
        !is.finite(pso_info$early_stopping) ||
        pso_info$early_stopping < 1 ||
        pso_info$early_stopping != as.integer(pso_info$early_stopping)) {
        stop("`early_stopping` must be a positive integer.")
    }
    
    if (pso_info$early_stopping > pso_info$max_iter) {
        stop("`early_stopping` cannot exceed `max_iter`.")
    }
    
    # Convergence tolerance
    if (!is.numeric(pso_info$tol) ||
        length(pso_info$tol) != 1 ||
        !is.finite(pso_info$tol) ||
        pso_info$tol <= 0) {
        stop("`tol` must be a positive numeric value.")
    }
    
    # Cognitive coefficient
    if (!is.numeric(pso_info$c1) ||
        length(pso_info$c1) != 1 ||
        !is.finite(pso_info$c1) ||
        pso_info$c1 < 0) {
        stop("`c1` must be a non-negative numeric value.")
    }
    
    # Social coefficient
    if (!is.numeric(pso_info$c2) ||
        length(pso_info$c2) != 1 ||
        !is.finite(pso_info$c2) ||
        pso_info$c2 < 0) {
        stop("`c2` must be a non-negative numeric value.")
    }
    
    # Starting inertia weight
    if (!is.numeric(pso_info$w0) ||
        length(pso_info$w0) != 1 ||
        !is.finite(pso_info$w0) ||
        pso_info$w0 <= 0) {
        stop("`w0` must be a positive numeric value.")
    }
    
    # Ending inertia weight
    if (!is.numeric(pso_info$w1) ||
        length(pso_info$w1) != 1 ||
        !is.finite(pso_info$w1) ||
        pso_info$w1 <= 0) {
        stop("`w1` must be a positive numeric value.")
    }
    
    if (pso_info$w1 > pso_info$w0) {
        stop("`w1` must not exceed `pso_info$w0`.")
    }
    
    # Inertia weight variation
    if (!is.numeric(pso_info$w_var) ||
        length(pso_info$w_var) != 1 ||
        !is.finite(pso_info$w_var) ||
        pso_info$w_var < 0 ||
        pso_info$w_var > 1) {
        stop("`w_var` must be between 0 and 1.")
    }
    
    # Velocity clamping factor
    if (!is.numeric(pso_info$vk) ||
        length(pso_info$vk) != 1 ||
        !is.finite(pso_info$vk) ||
        pso_info$vk <= 0) {
        stop("`vk` must be a positive numeric value.")
    }
    
    invisible(TRUE)

}


check_init_values <- function(init_values) {
    
    if (!inherits(init_values, "InitialValue")) {
        stop("`init_values` must be an object of class `InitialValue`.")
    }
    
    # Initial particle positions
    if (!is.null(init_values$init_swarm)) {
        
        if (!is.matrix(init_values$init_swarm)) {
            stop("`init_swarm` must be a matrix.")
        }
        
        if (any(!is.finite(init_values$init_swarm))) {
            stop("`init_swarm` must contain only finite values.")
        }
    }
    
    # Initial locally optimal design
    if (!is.null(init_values$init_local)) {
        
        if (!is.numeric(init_values$init_local) || !is.vector(init_values$init_local)) {
            stop("`init_local` must be a numeric vector.")
        }
        
        if (any(!is.finite(init_values$init_local))) {
            stop("`init_local` must contain only finite values.")
        }
        
        if (length(init_values$init_local) != 3) {
            stop("`init_local` must have length 3.")
        }
    }
    
    
    # Initial coefficient matrix
    if (!is.null(init_values$init_coef_mat)) {
        
        if (!is.matrix(init_values$init_coef_mat)) {
            stop("`init_coef_mat` must be a matrix.")
        }
        
        if (nrow(init_values$init_coef_mat) < 1) {
            stop("`init_coef_mat` must have at least one row.")
        }
        
        if (any(!is.finite(init_values$init_coef_mat))) {
            stop("`init_coef_mat` must contain only finite values.")
        }
    }
    
    invisible(TRUE)
    
}


check_minimax_alt_args <- function(
        design_type,
        distribution,
        design_info,
        pso_info,
        coef,
        coef_lower,
        coef_upper,
        init_values,
        highest_level,
        n_threads,
        verbose,
        seed
) {
    
    # init_values
    if (!is.null(init_values) &&
        !inherits(init_values, "InitialValue")) {
        stop("`init_values` must be an `InitialValue` object.")
    }
    
    
    # design_info
    if (!inherits(design_info, "DesignInfo")) {
        stop("`design_info` must be a `DesignInfo` object.")
    }
    
    
    # pso_info
    if (!inherits(pso_info, "PSOInfo")) {
        stop("`pso_info` must be a `PSOInfo` object.")
    }
    
    
    # design_type
    if (length(design_type) != 1 || !design_type %in% c("locally", "minimax")) {
        stop("`design_type` must be either `locally` or `minimax`.")
    }
    
    
    # distribution
    if (length(distribution) != 1 || 
        !distribution %in% c("weibull", "lognormal")) {
        stop("`distribution` must be `weibull` or `lognormal`.")
    }
    
    # coef / coef_lower / coef_upper
    if (design_type == "locally") {
        
        if (is.null(coef)) {
            stop("`coef` must be provided when `design_type = locally`.")
        }
        
        if (!is.numeric(coef) || any(!is.finite(coef))) {
            stop("`coef` must be a numeric vector containing only finite values.")
        }
        
        if (length(coef) - 1 != design_info$n_factor) {
            stop(
                "`coef` must contain 1 intercept and ",
                design_info$n_factor,
                " factor coefficients, ",
                "for a total of ",
                design_info$n_factor + 1,
                " elements."
            )
        }
        
    } else {
        
        if (is.null(coef_lower)) {
            stop("`coef_lower` must be provided when `design_type = minimax`.")
        }
        
        if (is.null(coef_upper)) {
            stop("`coef_upper` must be provided when `design_type = minimax`.")
        }
        
        if (!is.numeric(coef_lower) || any(!is.finite(coef_lower))) {
            stop(
                "`coef_lower` must be a numeric vector containing only finite values."
            )
        }
        
        if (!is.numeric(coef_upper) || any(!is.finite(coef_upper))) {
            stop(
                "`coef_upper` must be a numeric vector containing only finite values."
            )
        }
        
        if (length(coef_lower) - 1 != design_info$n_factor) {
            stop(
                "`coef_lower` must contain 1 intercept and ",
                design_info$n_factor,
                " factor coefficients, ",
                "for a total of ",
                design_info$n_factor + 1,
                " elements."
            )
        }
        
        if (length(coef_lower) != length(coef_upper)) {
            stop("`coef_lower` and `coef_upper` must have the same length.")
        }
        
        if (any(coef_lower > coef_upper)) {
            stop("Every element of `coef_lower` must be smaller than `coef_upper`.")
        }
        
        if(!is.null(coef)) { 
            
            if (!is.numeric(coef) || any(!is.finite(coef))) {
                stop("`coef` must be a numeric vector containing only finite values.")
            }
            
            if (length(coef) - 1 != design_info$n_factor) {
                stop(
                    "`coef` must contain 1 intercept and ",
                    design_info$n_factor,
                    " factor coefficients, ",
                    "for a total of ",
                    design_info$n_factor + 1,
                    " elements."
                )
            }
            
            if (any(coef < coef_lower) || any(coef > coef_upper)) {
                stop("Every element of `coef` must be within its corresponding bounds.")
            }
            
        }
    
    }

    
    # highest_level
    if (!is.logical(highest_level) ||
        length(highest_level) != 1 ||
        is.na(highest_level)) {
        stop("`highest_level` must be a single TRUE or FALSE value.")
    }
    
    
    # n_threads
    if (!is.numeric(n_threads) ||
        length(n_threads) != 1 ||
        !is.finite(n_threads) ||
        n_threads < 1 ||
        n_threads != as.integer(n_threads)) {
        stop("`n_threads` must be a positive integer.")
    }
    
    
    # verbose
    if (!is.logical(verbose) ||
        length(verbose) != 1 ||
        is.na(verbose)) {
        stop("`verbose` must be a single TRUE or FALSE value.")
    }
    
    
    # seed
    if (!is.numeric(seed) ||
        length(seed) != 1 ||
        !is.finite(seed) ||
        seed != as.integer(seed)) {
        stop("`seed` must be an integer.")
    }
    
    
    # relationship between reparam and design_type
    if (!design_info$reparam && design_type != "locally") {
        stop(
            "Non-reparameterization is only available for `design_type = locally`."
        )
    }
    
    # init_values
    if (!is.null(init_values)) {
        
        if (!is.null(init_values$init_swarm)) {
            
            d_swarm <- design_info$n_support * design_info$n_factor + 
                design_info$n_support - 1
            
            var_lower <- c(rep(design_info$x_l, design_info$n_support * design_info$n_factor), 
                           rep(0, design_info$n_support - 1))
            
            var_upper <- c(rep(design_info$x_h, design_info$n_support * design_info$n_factor), 
                           rep(1, design_info$n_support - 1))
            
            if (nrow(init_values$init_swarm) != pso_info$n_swarm) {
                stop(
                    "`init_swarm` must have ", pso_info$n_swarm,
                    " rows (one row for each particle)."
                )
            }
            
            if (ncol(init_values$init_swarm) != d_swarm) {
                stop(
                    "`init_swarm` must have ", d_swarm,
                    " columns (number of stress levels times factors and their proportions)."
                )
            }
            
            for (i in 1:nrow(init_values$init_swarm)) {
                
                if (any(var_upper < init_values$init_swarm[i,]) || 
                    any(init_values$init_swarm[i,] < var_lower)) {
                    
                    print(var_upper)
                    print(var_lower)
                    stop("Initial particle ", i,
                        " contains values outside the specified bounds.")
                }
            }
        }
        
        if (!is.null(init_values$init_local)) {
            
            local_lower <- c(rep(design_info$x_l, 2), 0)
            local_upper <- c(rep(design_info$x_h, 2), 1)
            
            if (any(init_values$init_local < local_lower) ||
                any(init_values$init_local > local_upper)) {
                    
                stop("Every element of `init_local` must be within its corresponding bounds.")
            }
        }
        
        if (!is.null(init_values$init_coef_mat)) {
            
            if (ncol(init_values$init_coef_mat) != design_info$n_factor + 1) {
                stop(
                    "`init_coef_mat` must have ", design_info$n_factor + 1,
                    " columns."
                )
            }
        }
        
    }
    
    invisible(TRUE)
}


check_design_measure <- function(best_design, design_info) {
    
    if (!inherits(design_info, "DesignInfo")) {
        stop("`design_info` must be a `DesignInfo` object.")
    }
    
    if (!is.matrix(best_design)) {
        stop("`best_design` must be a matrix.")
    }
    
    if (design_info$n_factor != nrow(best_design) - 1) {
        stop(
            "`best_design` has ", nrow(best_design) - 1,
            " stress factors, but `j_factor` is ",
            design_info$n_factor, "."
        )
    }
    
    if (design_info$n_support != ncol(best_design)) {
        stop(
            "`best_design` has ", ncol(best_design),
            " levels, but `k_levels` is ",
            design_info$n_support, "."
        )
    }
    
    if (any(best_design[nrow(best_design), ] < 0) ||
        any(best_design[nrow(best_design), ] > 1)) {
        stop(
            "The proportions of `best_design` must be between 0 and 1."
        )
    }
    
    if (!isTRUE(all.equal(sum(best_design[nrow(best_design), ]), 1))) {
        stop(
            "The proportions of `best_design` must sum to 1."
        )
    }
    
    invisible(TRUE)
}


check_model_set <- function(model_set, design_info) {
    
    if (!inherits(design_info, "DesignInfo")) {
        stop("`design_info` must be a `DesignInfo` object.")
    }
    
    if (!is.matrix(model_set)) {
        stop("`model_sel` must be a matrix.")
    }
    
    if (nrow(model_set) < 1) {
        stop("`init_coef_mat` must have at least one row.")
    }
    
    
    # coefficients
    coef_mat <- matrix(
        as.numeric(model_set[, 1:(ncol(model_set) - 1)]),
        nrow = nrow(model_set[, 1:(ncol(model_set) - 1)])
    )
    
    if (any(!is.finite(coef_mat))) {
        stop("Coefficients of `model_set` must contain only finite values.")
    }

    if (ncol(coef_mat) != design_info$n_factor + 1) {
        stop(
            "Coefficients of `model_set` must have ", design_info$n_factor + 1,
            " columns."
        )
    }
    
    
    # distribution
    dist_vec <- model_set[, ncol(model_set)]
    
    if (!all(dist_vec %in% c("weibull", "lognormal"))) {
        stop("`distribution` must be `weibull` or `lognormal`.")
    }
    
    invisible(TRUE)
    
}

