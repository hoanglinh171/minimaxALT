#' @export
summary.OptimalALT <- function(object, ...) {
    stopifnot(inherits(object, "OptimalALT"))
    
    design <- extract_design(object)
    
    cat("Summary of generated optimal ALT design\n")
    cat("-----------------------------------------------\n")
    cat("X: Stress levels, W: Corresponding proportion\n")
    print(design)
    cat("\nObjective Value:", object$fg_best, "\n")
    cat("Max directional derivative:", object$max_directional_derivative, "\n")
    cat("\n")
    
    invisible(object)
}

