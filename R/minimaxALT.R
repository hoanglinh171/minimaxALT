#' Find Optimal ALT Design Using Hybrid Algorithm
#'
#' Runs hybrid algorithm combining PSO and Nelder-Mead to find the optimal design of accelerated life test (ALT).
#'
#' @param design_type Character. One of c("locally", "minimax").
#' @param distribution Character. Failure distribution, one of 
#'  c("weibull, "lognormal").
#' @param design_info A `DesignInfo` object from \code{\link{set_design_info}} containing 
#'  design specifications.
#' @param pso_info A `PSOInfo` object from \code{\link{pso_setting}} defining 
#'  PSO hyperparameters.
#' @param coef Optional. Fixed model coefficients. 
#'  Required if \code{design_type = "locally"}.
#' @param coef_lower Optional. Lower bounds for model parameters. 
#'  Required if \code{design_type = "minimax"}.
#' @param coef_upper Optional. Upper bounds for model parameters. 
#'  Required if \code{design_type = "minimax"}.
#' @param init_values Optional. An `InitialValue` object of initial values 
#'  from \code{\link{initialize_values}}.
#' @param highest_level Logical. Whether the highest stress level of the 
#'  generated design is the upper bound of stress range \code{x_h}. 
#'  Default value is \code{FALSE}.
#' @param n_threads Integer. Number of threads for parallel processing.
#' @param verbose Logical. If \code{TRUE}, print optimization progress.
#' @param seed Integer. Seed for reproducibility
#' 
#' @return An object of class OptimalALT
#' 
#' @examples
#' 
#' design_info <- set_design_info(
#'     k_levels=3, 
#'     j_factor=1, 
#'     n_unit=300, 
#'     censor_time=183, 
#'     p=0.1, 
#'     use_cond=0, 
#'     sigma=0.6
#' )
#' 
#' pso_info <- pso_setting(
#'     n_swarm=5, 
#'     max_iter=10, 
#'     early_stopping=5, 
#'     tol=0.00001
#' )
#' 
#' res <- find_optimal_alt(
#'     design_type="minimax", 
#'     distribution="weibull", 
#'     design_info=design_info, 
#'     pso_info=pso_info, 
#'     coef_lower=c(10^-6, 0.7),
#'     coef_upper=c(10^-3, 0.99), 
#'     highest_level = TRUE,
#'     verbose = FALSE,
#'     n_threads = 1
#' )
#' 
#' print(res)
#' summary(res)
#' plot(res, x_l=0, x_h=1)
#' 
#' @references 
#' \enumerate{
#'   \item Chen P (2024). _globpso: Particle Swarm Optimization Algorithms and Differential Evolution for Minimization Problems_. R package version 1.2.1, <https://github.com/PingYangChen/globpso>.
#'   \item Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization. In Proceedings of the IEEE International Conference on Neural Networks (ICNN) (Vol. 4, pp. 1942–1948).
#'   \item Lee, I. C., Chen, R. B., Wong, W. K., (in press). Optimal Robust Strategies for Accelerated Life Tests and Fatigue Testing of Polymer Composite Materials. Annals of Applied Statistics. <https://imstat.org/journals-and-publications/annals-of-applied-statistics/annals-of-applied-statistics-next-issues/>
#'   \item Meeker, W. Q., & Escobar, L. A. (1998). Statistical methods for reliability data. New York: Wiley-Interscience.
#'   \item Nelder, J. A. and Mead, R. (1965). A simplex algorithm for function minimization. Computer Journal, 7, 308--313. 10.1093/comjnl/7.4.308.
#' }
#' 
#' @export
find_optimal_alt <- function(design_type, distribution,
                             design_info, pso_info,
                             coef = NULL,
                             coef_lower = NULL, coef_upper = NULL,
                             init_values = NULL,
                             highest_level = FALSE,
                             n_threads = 1,
                             verbose = TRUE,
                             seed = 42) {
  
    ## condition check
    args <- mget(names(formals()))
    do.call(check_minimax_alt_args, args)
    
    
    ## number of threads
    max_cores <- parallel::detectCores()
    if (n_threads >= max_cores) {
        
        warning(sprintf("Number of threads available is %d, ", max_cores),
                sprintf("Using %d threads instead.", max_cores))
        
        n_threads <-  round(max_cores, digits = 0)
    }
    
    
    ## seed
    seed <- round(seed, digits = 0)
    set.seed(seed)
    
    
    ## coef pseudo-bounds for locally optimal design
    if (design_type == "locally") {

        coef_upper <- rep(10000, design_info$n_factor + 1)
        coef_lower <- rep(-10000, design_info$n_factor + 1)
    
    }
    
    
    ## bounds of swarm
    x_l <- design_info$x_l
    x_h <- design_info$x_h
    n_support <- design_info$n_support
    n_factor <- design_info$n_factor
    
    var_lower <- c(rep(x_l, n_support * n_factor), rep(0, n_support - 1))
    var_upper <- c(rep(x_h, n_support * n_factor), rep(1, n_support - 1))
  
    if (highest_level) {
        var_lower[seq(1, n_support * n_factor, by = n_support)] <- x_h
    }

    
    ## particle size
    d_swarm <- length(var_lower)
    n_swarm <- pso_info$n_swarm
    
    
    ## bounds of nelder-mead
    local_lower <- c(rep(design_info$x_l, 2), 0)
    local_upper <- c(rep(design_info$x_h, 2), 1)
  
    
    ## get initial values
    init_coef <- NULL
    init_swarm <- NULL
    init_local <- NULL
    init_coef_mat <- NULL
    
    if(!is.null(coef)) {
        init_coef <- coef
    } else {
        init_coef <- runif(n_factor + 1, min = -40, max = 40)
    }
    
    if (!is.null(init_values)) {
        
        init_swarm <- init_values$init_swarms
        init_local <- init_values$init_local
        init_coef_mat <- init_values$init_coef_mat
    
    } 
      
    if (is.null(init_swarm)) {
        
        init_swarm <- matrix(runif(d_swarm * n_swarm), ncol = d_swarm)
        init_swarm <- init_swarm * 
            matrix(rep(var_upper - var_lower, n_swarm), ncol = d_swarm, byrow = TRUE) + 
            matrix(rep(var_lower, n_swarm), ncol = d_swarm, byrow = TRUE)
    
    } 

    if(is.null(init_local)) {
        init_local <- c(1, 0.6, 0.3)
    }
    
    if(is.null(init_coef_mat)) {
        
        init_coef_mat <- 10 * as.matrix(expand.grid(rep(list(c(1, -1)), n_factor + 1)))
        # init_coef_mat = rbind(init_coef_mat, 40 * as.matrix(expand.grid(rep(list(c(1, -1)), n_factor + 1))))
        
    } else {
        
        for(i in 1:nrow(init_coef_mat)) {
            init_coef_mat[i, ] <- get_outbound_sigmoid(init_coef_mat[i, ], coef_lower, coef_upper)
        }
    }
    
    ## define pso info
    pso_info$var_upper <- var_upper
    pso_info$var_lower <- var_lower
    pso_info$d_swarm <- d_swarm
    pso_info$init_swarm <- t(init_swarm)

    ## inner parameters
    init_bound_info <- list()
    init_bound_info$opt_coef <- init_coef
    init_bound_info$coef_upper <- coef_upper
    init_bound_info$coef_lower <- coef_lower
    init_bound_info$init_coef <- init_coef
    init_bound_info$opt_local <- init_local
    init_bound_info$local_upper <- local_upper
    init_bound_info$local_lower <- local_lower
    init_bound_info$init_local <- init_local
    init_bound_info$model <- ifelse(distribution == "weibull", 1, 
                                    ifelse(distribution == "lognormal", 2, 3))
    init_bound_info$opt_distribution <- init_bound_info$model
  
    # nelder-mean
    nelder_mead_settings = list()
    nelder_mead_settings$init_coef_mat = t(init_coef_mat)
    
    # run
    design_type_code <- ifelse(design_type == "locally", 1, 2)
    minimax_design <- minimax_alt(design_type_code, pso_info, design_info, init_bound_info,
                                nelder_mead_settings,
                                n_threads, verbose, seed)

  
    class(minimax_design) <- "OptimalALT"
  
    return(minimax_design)

}

