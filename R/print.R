#' Print an OptimalALT Object 
#' 
#' Prints a summary of the optimization results and optimality check 
#'  for an \code{OptimalALT} object. 
#' 
#' @param x An object of class \code{OptimalALT}, typically returned 
#'  by \code{\link{find_optimal_alt}}. 
#' @param ... Additional arguments.
#'  
#' @return Invisibly returns the original \code{OptimalALT} object. 
#' 
#' @examples 
#' \dontrun{ 
#' result <- find_optimal_alt(...) 
#' print(result) 
#' }
#' 
#' @export
print.OptimalALT <- function(x, ...) {
    stopifnot(inherits(x, "OptimalALT"))
    
    iterations <- length(x$fg_best_hist) - 1
    model_num <- nrow(x$model_set)
    
    cat("PSO implementation\n")
    cat("-----------------------------------------------\n")
    cat("Iterations:", iterations, "\n")
    cat("Optimal particle:", x$g_best, "\n")
    cat("xive value:", x$fg_best, "\n")
    cat("\n")
    
    cat("Optimality check\n")
    cat("-----------------------------------------------\n")
    cat("Number of model candidates:", model_num,"\n")
    cat("Max directional derivative:", x$max_directional_derivative, "\n")
    cat("\n")
    
    invisible(x)
}


#' Print a DesignInfo Object
#' 
#' Prints the design specifications stored in a \code{DesignInfo} object.
#' 
#' @param x An object of class \code{DesignInfo}, typically returned
#'  by \code{\link{set_design_info}}.
#' @param ... Additional arguments.
#'  
#' @return Invisibly returns the original \code{DesignInfo} object. 
#' 
#' @examples
#' \dontrun{
#' design_info <- set_design_info(...)
#' print(design_info)
#' } 
#' 
#' @export
print.DesignInfo <- function(x, ...) {
    stopifnot(inherits(x, "DesignInfo"))
    
    cat("ALT design specifications\n")
    cat("-----------------------------------------------\n")
    cat("Optimality type: C-optimality\n")
    cat("Number of stress levels:", x$n_support, "\n")
    cat("Number of factors:", x$n_factor, "\n")
    cat("Number of testing units:", x$n_unit, "\n")
    cat("Censoring time:", x$censor_time, "\n")
    cat("Stress level at the use condition:", 
        paste0("(", paste(x$use_cond, collapse = ", "), ")"), "\n")
    cat("Stress range: [", x$x_l, ", ", x$x_h, "]\n", sep = "")
    cat("Lifetime percentile to be estimated at the use condition:", x$p, "\n")
    cat("Coefficients are reparameterized into failure probability?", x$reparam, "\n")
    cat("Scale parameter:", x$sigma, "\n")
    cat("\n")
    
    invisible(x)
}


#' Print a PSOInfo Object 
#' 
#' Prints the particle swarm optimization (PSO) hyperparameters stored 
#' in a \code{PSOInfo} object.
#' 
#' @param x An object of class \code{PSOInfo}, typically returned 
#'  by \code{\link{pso_setting}}.
#' @param ... Additional arguments. 
#'  
#' @return Invisibly returns the original \code{PSOInfo} object. 
#'  
#' @examples
#' \dontrun{
#' pso_info <- pso_setting(...)
#' print(pso_info)
#' }
#'  
#' @export
print.PSOInfo <- function(x, ...) {
    stopifnot(inherits(x, "PSOInfo"))
    
    cat("PSO Hyperparameters\n")
    cat("-----------------------------------------------\n")
    cat("Number of particles:", x$n_swarm, "\n")
    cat("Maximum number of iterations:", x$max_iter, "\n")
    cat("Frequency of checking optimality:", x$early_stopping, "iterations\n")
    cat("Convergence tolerance:", x$tol, "\n")
    cat("Cognitive acceleration coefficient c1:", x$c1, "\n")
    cat("Social acceleration coefficient c2:", x$c2, "\n")
    cat("Starting inertia weight w0:", x$w0, "\n")
    cat("Ending inertia weight w1:", x$w1, "\n")
    cat("Inertia weight decay fraction:", x$w_var, "\n")
    cat("Velocity clamping factor:", x$vk, "\n")
    cat("\n")
    
    invisible(x)
    
}


#' Print an InitialValue Object
#' 
#' Prints a concise summary of the initial values used for the PSO and 
#' Nelder-Mead optimization procedures. 
#' 
#' @param x An object of class \code{InitialValue}, typically created 
#'  by \code{\link{initialize_values}}.
#' @param ... Additional arguments.
#' 
#' @return Invisibly returns the \code{InitialValue} object. 
#' 
#' @examples
#' \dontrun{
#' init_values <- initialize_values(
#'     init_swarm = init_swarm,
#'     init_local = init_local,
#'     init_coef = init_coef
#' )
#' 
#' print(init_values)
#' }
#' 
#' @export
print.InitialValue <- function(x, ...) {
    stopifnot(inherits(x, "InitialValue"))
    
    cat("Initialized Values for PSO and Nelder-Mead Optimization\n")
    cat("-----------------------------------------------\n")
    cat("Number of initial particles for PSO:", nrow(x$init_swarm), "\n")
    cat("Particle size of PSO:", ncol(x$init_swarm), "\n")
    cat("Initial value (design) for Nelder-Mead optimization:", 
        paste0("(", paste(x$init_local, collapse = ", "), ")"), "\n")
    cat("Multi-start Nelder-Mead optimization:", nrow(x$init_coef), "times\n")
    cat("\n")
    
    invisible(x)
}


#' Print an OptimalityCheck Object
#' 
#' Prints the results of an optimality check for an accelerated life
#'  test (ALT) design.
#'  
#' @param x An object of class \code{OptimalityCheck}, typically 
#'  returned by \code{\link{check_optimality}}. 
#' @param show_candidates Logical. If \code{TRUE}, the candidate model 
#'  set used in the optimality check is printed. Default is \code{FALSE}.
#' @param ... Additional arguments.
#'  
#' @return Invisibly returns the original \code{OptimalityCheck} object. 
#' 
#' @examples
#' \dontrun{
#' optimality <- check_optimality(...)
#' print(optimality)
#' }
#' 
#' @export
print.OptimalityCheck <- function(x, show_candidates = FALSE, ...) {
    stopifnot(inherits(x, "OptimalityCheck"))
    
    model_num <- nrow(x$model_set)
    colnames(x$model_set) <- c("coef_1", "coef_2", "distribution")
    
    cat("Optimality check\n")
    cat("-----------------------------------------------\n")
    cat("Design: \n")
    print(x$design)
    cat("\n")
    cat("Number of model candidates:", model_num,"\n")
    if(show_candidates) {
        print(x$model_set)
        cat("\n")
    }
    cat("Max directional derivative:", x$max_directional_derivative, "\n")
    cat("\n")
    
    invisible(x)
}



