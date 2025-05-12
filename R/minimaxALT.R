

## Find optimal design ----------------------------------------------------------
find_optimal_alt <- function(design_type, distribution,
                             design_info, pso_info,
                             coef = NULL,
                             coef_lower = NULL, coef_upper = NULL,
                             init_values = NULL,
                             highest_level = TRUE,
                             n_threads = round(detectCores() * 0.5, digits = 0),
                             verbose = TRUE) {
  
  
  ## Condition check
  if (design_type == 1) {
    stopifnot(!is.null(coef), design_info$n_factor == length(coef) - 1)
    coef_upper = rep(10000, design_info$n_factor + 1)
    coef_lower = rep(-10000, design_info$n_factor + 1)
  } else if (design_type == 2) {
    stopifnot(!is.null(coef_lower), !is.null(coef_upper), 
              length(coef_upper) == length(coef_lower), 
              design_info$n_factor == length(coef_lower) - 1)
  } else {
    stop("Design type must be 1 (locally optimal design) or 2 (minimax design)")
  }
  
  stopifnot(distribution == 1 || distribution == 2 || distribution == 3)
  
  stopifnot(design_info$opt_type == "C" || design_info$opt_type == "D")
  if (design_info$opt_type == "D") {
    stop("D-optimal design is under development.")
  }
  
  stopifnot(is.numeric(design_info$n_support), is.numeric(design_info$n_factor), 
            is.numeric(design_info$n_unit), 
            is.numeric(design_info$censor_time), is.numeric(design_info$sigma), 
            is.numeric(design_info$p), 
            is.numeric(design_info$x_l), is.numeric(design_info$x_h))
  
  stopifnot(is.logical(design_info$reparam), is.logical(design_info$degenerate), 
            is.logical(verbose))
  
  use_cond = c(design_info$use_cond)
  stopifnot(design_info$n_factor == length(use_cond))
  
  max_cores = parallel::detectCores()
  if (n_threads > 0.8 * max_cores) {
    cat("Number of threads available is", max_cores, ".\n")
    cat("It is recommended to run with ", round(0.8 * max_cores, digits = 0), "threads at most.\n")
    n_threads = round(0.8 * max_cores, digits = 0)
  }
  
  
  ## Define design info
  design_info$use_cond = use_cond
  
  ## Define bounds of swarm
  x_l = design_info$x_l
  x_h = design_info$x_h
  n_support = design_info$n_support
  n_factor = design_info$n_factor
  
  
  var_lower <- c(rep(x_l, n_support * n_factor), rep(0, n_support - 1))
  var_upper <- c(rep(x_h, n_support * n_factor), rep(1, n_support - 1))
  
  if (highest_level) {
    var_lower[seq(1, n_support * n_factor, by = n_support)] <- x_h
  }

    
  d_swarm = length(var_lower)
  
  
  ## Get pso info
  stopifnot(all(names(pso_info) == names(pso_setting())))
  n_swarm <- pso_info$n_swarm
  
  
  ## Get initial values
  init_swarm = NULL
  init_coef = NULL
  init_local = NULL
  init_coef_mat = NULL
  
  if (!is.null(init_values)) {

    stopifnot(all(names(init_values) == names(initialize_values())))
    
    init_swarm = init_values$init_swarms
    
    if(!is.null(coef)) {
      init_coef = coef
    }
    
    init_local = init_values$init_local
    init_coef_mat = init_values$init_coef_mat
  } 
  
  
  if (is.null(init_swarm)) {
    init_swarm = matrix(runif(d_swarm * n_swarm), ncol = d_swarm)
    init_swarm = init_swarm * matrix(rep(var_upper - var_lower, n_swarm), ncol = d_swarm, byrow = TRUE) + 
      matrix(rep(var_lower, n_swarm), ncol = d_swarm, byrow = TRUE)
  } else {
    stopifnot(is.matrix(init_swarm),
              ncol(init_swarm) == d_swarm,
              nrow(init_swarm) == n_swarm,
              all(is.finite(init_swarm))
    )
    
    for (i in 1:nrow(init_swarm)) {
      stopifnot(all(var_upper >= init_swarm[i,]), 
                all(init_swarm[i,] >= var_lower)
      )
    }
  }
  
  
  if(is.null(init_coef)) {
    if(is.null(coef)) {
      init_coef = runif(n_factor + 1, min = -40, max = 40)
    } else {
      init_coef = coef
    }

    
  } else {
    stopifnot(length(init_coef) == n_factor + 1, 
              all(coef_upper >= init_coef), 
              all(is.finite(init_coef)), 
              all(init_coef >= coef_lower)
    )
  }
  
  
  local_lower = c(rep(x_l, 2), 0)
  local_upper = c(rep(x_h, 2), 1)
  
  if(is.null(init_local)) {
    init_local = c(1, 0.6, 0.3)
  } else {
    stopifnot(length(local_lower) == length(init_local), 
              all(local_upper >= init_local), 
              all(is.finite(init_local)), 
              all(init_local >= local_lower)
    )
  }
  
  
  
  
  if(is.null(init_coef_mat)) {
    init_coef_mat = 10 * as.matrix(expand.grid(rep(list(c(1, -1)), n_factor + 1)))
    # init_coef_mat = rbind(init_coef_mat, 40 * as.matrix(expand.grid(rep(list(c(1, -1)), n_factor + 1))))
  } else {
    stopifnot(is.matrix(init_coef_mat),
              ncol(init_coef_mat) == n_factor + 1,
              all(is.finite(init_coef_mat))
    )
    
    for(i in 1:nrow(init_coef_mat)) {
      init_coef_mat[i, ] = get_outbound_sigmoid(init_coef_mat[i, ], coef_lower, coef_upper)
    }
  }
    
  
  
  
  ## Define pso info
  pso_info$var_upper <- var_upper
  pso_info$var_lower <- var_lower
  pso_info$d_swarm <- d_swarm
  pso_info$init_swarm <- t(init_swarm)

  
  ## Inner parameters
  init_bound_info <- list()
  init_bound_info$opt_coef = init_coef
  init_bound_info$coef_upper = coef_upper
  init_bound_info$coef_lower = coef_lower
  init_bound_info$init_coef = init_coef
  init_bound_info$opt_local = init_local
  init_bound_info$local_upper = local_upper
  init_bound_info$local_lower = local_lower
  init_bound_info$init_local = init_local
  init_bound_info$model = distribution
  init_bound_info$opt_distribution = distribution
  
  
  nelder_mead_settings = list()
  nelder_mead_settings$init_coef_mat = t(init_coef_mat)
  

  minimax_design <- minimax_alt(design_type, pso_info, design_info, init_bound_info,
                                nelder_mead_settings,
                                n_threads, verbose)

  
  class(minimax_design) <- "OptimalALT"
  return(minimax_design)
}

