#' Extract Optimal Design
#' 
#' Extract design from an OptimalALT object
#'  
#' @param optimal_alt An object of class `OptimalALT`.
#' 
#' @return A matrix containing stress levels and allocated proportion of 
#'  the design, each row corresponds to stress levels of each factor, 
#'  while the last row is their proportions.
#' 
#' @examples 
#' 
#' design_info <- set_design_info(
#'     k_levels=2, 
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
#'     tol=0.0001
#' )
#' 
#' example1_locally <- find_optimal_alt(
#'     design_type="locally", 
#'     distribution="weibull",
#'     design_info=design_info, 
#'     pso_info=pso_info, 
#'     coef = c(0.001, 0.9),
#'     highest_level=TRUE,
#'     verbose=TRUE, 
#'     n_threads = 1)
#'     
#' opt_design <- extract_design(example1_locally)   
#' 
#'   
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
