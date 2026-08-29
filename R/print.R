#' @export
print.OptimalALT <- function(object, ...) {
    stopifnot(inherits(object, "OptimalALT"))
    
    iterations <- length(object$fg_best_hist) - 1
    model_num <- nrow(object$model_set)
    
    cat("PSO implementation\n")
    cat("-----------------------------------------------\n")
    cat("Iterations:", iterations, "\n")
    cat("Optimal particle:", object$g_best, "\n")
    cat("Objective value:", object$fg_best, "\n")
    cat("\n")
    
    cat("Optimality check\n")
    cat("-----------------------------------------------\n")
    cat("Number of model candidates:", model_num,"\n")
    cat("Max directional derivative:", object$max_directional_derivative, "\n")
    cat("\n")
    
    invisible(object)
}


#' @export
print.DesignInfo <- function(object) {
    stopifnot(inherits(object, "DesignInfo"))
    
    cat("ALT design specifications\n")
    cat("-----------------------------------------------\n")
    cat("Number of stress levels:", object$n_support, "\n")
    cat("Number of factors:", object$n_factor, "\n")
    cat("Number of testing units:", object$n_unit, "\n")
    cat("Censoring time:", object$censor_time, "\n")
    cat("Stress level at the use condition:", 
        paste0("(", paste(object$use_cond, collapse = ", "), ")"), "\n")
    cat("Stress range: [", object$x_l, ", ", object$x_h, "]\n", sep = "")
    cat("Lifetime percentile to be estimated at the use condition:", object$p, "\n")
    cat("Coefficients are reparameterized into failure probability?", object$reparam, "\n")
    cat("Scale parameter:", object$sigma, "\n")
    cat("\n")
    
    invisible(object)
}


#' @export
print.PSOInfo <- function(object) {
    stopifnot(inherits(object, "PSOInfo"))
    
    cat("PSO Hyperparameters\n")
    cat("-----------------------------------------------\n")
    cat("Number of particles:", object$n_swarm, "\n")
    cat("Maximum number of iterations:", object$max_iter, "\n")
    cat("Frequency of checking optimality:", object$early_stopping, "iterations\n")
    cat("Convergence tolerance:", object$tol, "\n")
    cat("Cognitive acceleration coefficient c1:", object$c1, "\n")
    cat("Social acceleration coefficient c2:", object$c2, "\n")
    cat("Starting inertia weight w0:", object$w0, "\n")
    cat("Ending inertia weight w1:", object$w1, "\n")
    cat("Inertia weight decay fraction:", object$w_var, "\n")
    cat("Velocity clamping factor:", object$vk, "\n")
    cat("\n")
    
    invisible(object)
    
}


#' @export
print.InitialValue <- function(object) {
    stopifnot(inherits(object, "InitialValue"))
    
    cat("Initialized Values for PSO and Nelder-Mead Optimization\n")
    cat("-----------------------------------------------\n")
    cat("Number of initial particles for PSO:", nrow(object$init_swarm), "\n")
    cat("Particle size of PSO:", ncol(object$init_swarm), "\n")
    cat("Initial value (design) for Nelder-Mead optimization:", 
        paste0("(", paste(object$init_local, collapse = ", "), ")"), "\n")
    cat("Multi-start Nelder-Mead optimization:", nrow(object$init_coef), "times\n")
    cat("\n")
    
    invisible(object)
}


#' @export
print.OptimalityCheck <- function(object, show_candidates = FALSE) {
    stopifnot(inherits(object, "OptimalityCheck"))
    
    model_num <- nrow(object$model_set)
    colnames(object$model_set) <- c("coef_1", "coef_2", "distribution")
    
    cat("Optimality check\n")
    cat("-----------------------------------------------\n")
    cat("Design: \n")
    print(object$design)
    cat("\n")
    cat("Number of model candidates:", model_num,"\n")
    if(show_candidates) {
        print(object$model_set)
        cat("\n")
    }
    cat("Max directional derivative:", object$max_directional_derivative, "\n")
    cat("\n")
    
    invisible(object)
}



