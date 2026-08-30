#' Summarize an OptimalALT Object
#' 
#' Provides a concise summary of the generated optimal accelerated life
#' test (ALT) design and its optimality-check results. 
#' 
#' @param object An object of class \code{OptimalALT}, typically returned 
#'  by \code{\link{find_optimal_alt}}. 
#' @param ... Additional arguments.
#'  
#' @return Invisibly returns the original \code{OptimalALT} object.
#' 
#' @examples 
#' \dontrun{
#' result <- find_optimal_alt(...)
#' summary(result)
#' }
#' 
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

