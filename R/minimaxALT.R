#' Find Optimal ALT Design Using Hybrid Algorithm
#'
#' Runs hybrid algorithm combining PSO and Nelder-Mead to find the optimal design of accelerated life test (ALT).
#'
#' @param design_type Integer. 1: Locally optimal design, 2: Minimax design.
#' @param distribution Integer. The assumed failure time distribution, 1: Weibull, 2: Log-normal, 3: Model robust (both distribution Weibull and Log-normal).
#' @param design_info A DesignInfo object from `set_design_info()` containing design specifications.
#' @param pso_info A PSOInfo object from `pso_setting()` defining PSO hyperparameters.
#' @param coef Optional. Fixed model coefficients. Required if \code{design_type = 1}.
#' @param coef_lower Optional. Lower bounds for model parameters. Required if \code{design_type = 2}.
#' @param coef_upper Optional. Upper bounds for model parameters. Required if \code{design_type = 2}.
#' @param init_values Optional. An InitialValue object of initial values from `initialize_values()`.
#' @param highest_level Logical. Whether the highest stress level of the generated design is the upper bound of stress range \code{x_h}. Default value is \code{TRUE}.
#' @param n_threads Integer. Number of threads for parallel processing.
#' @param verbose Logical. If \code{TRUE}, print optimization progress.
#' @param seed Integer. Seed for reproducibility
#' @return
#' \describe{
#' \item{g_best}{The global best design found by the hybrid algorithm.}
#' \item{coef_best}{The parameters corresponding to the global best design.}
#' \item{distribution_best}{The distribution corresponding to the global best design.}
#' \item{max_directional_derivative}{Maximum directional derivative within design space, evaluated using equivalence theorem.}
#' \item{fg_best}{The objective function value corresponding to the global best design.}
#' \item{fg_best_hist}{A vector tracking the best objective function value of each iteration.}
#' \item{p_best}{A matrix containing each particle's personal best design found during the optimization.}
#' \item{fp_best}{A vector containing the objective function values corresponding to each particle's personal best.}
#' \item{g_hist}{All particle positions of each iteration.}
#' \item{coef_best_hist}{The parameters corresponding to the global best designs of each iteration.}
#' \item{distribution_best_hist}{The distribution corresponding to the global best designs of each iteration.}
#' \item{model_set}{A matrix containing distribution and model parameters of global best particles of each iteration, duplicated models are removed.}
#' \item{model_weight}{The weight assigned to each model in the model set.}
#' \item{equivalence_data}{Generated designs and their corresponding directional derivative given the optimal design \code{g_best}. Each design is a combination of factors with value in [0, 1]. These designs are data for plotting equivalence theorem plot.}
#' }
#' @examples
#' 
#' design_info <- set_design_info(k_levels=2, j_factor=1, n_unit=300, 
#'                                censor_time=183, p=0.1, use_cond=0, sigma=0.6)
#' 
#' pso_info <- pso_setting(n_swarm=32, max_iter=128, early_stopping=10, tol=0.01)
#' 
#' set.seed(42)
#' res <- find_optimal_alt(design_type=1, distribution=1, design_info=design_info, 
#'                         pso_info=pso_info, coef=c(0.001, 0.9), verbose = FALSE)
#' 
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
#' @name find_optimal_alt
#' @rdname find_optimal_alt
#' @importFrom Rcpp evalCpp cppFunction sourceCpp
#' @importFrom parallel detectCores
#' @importFrom stats runif
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

