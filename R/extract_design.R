#' Extract design from an OptimalALT object
#'
#' @param optimal_alt An object of class `OptimalALT`.
#' @return The extracted design.
#' @export
extract_design <- function(optimal_alt) {
    stopifnot(inherits(optimal_alt, "OptimalALT"))
    
    n_factor <- length(optimal_alt$coef_best) - 1
    n_support <- (length(optimal_alt$g_best) + 1) / (n_factor + 1)
    
    stress_levels <- matrix(optimal_alt$g_best[1:(n_support*n_factor)], 
                            ncol = n_support, byrow=TRUE)
    
    prop <- get_proportion(optimal_alt$g_best[(n_support*n_factor + 1):length(optimal_alt$g_best)])
    
    design <- rbind(stress_levels, prop)
    
    level_names <- paste0("X", 1:n_factor)
    level_names <- c(level_names, "W")
    rownames(design) <- level_names

    return(design)

}
