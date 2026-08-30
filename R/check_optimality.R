#' Check Optimality for Optimal Design
#'
#' Evaluates whether a design is optimal by the equivalence theorem.
#'
#' @param best_design A matrix containing stress levels and allocated 
#'  proportion of the design, each row corresponds to stress levels of 
#'  each factor, while the last row is their proportions.
#' @param model_set A matrix of models, including parameters and distribution, 
#'  that maximize the optimality criteria with the given best particle's position, 
#'  each columns corresponds to model coefficients, while the last column 
#'  is lifetime distribution.
#' @param design_info A `DesignInfo` object created by 
#'  function \code{\link{set_design_info}}.
#' @param seed Seed for reproducibility
#'
#' @return An OptimalityCheck object containing optimality check results
#' 
#' @examples
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
#' best_design <- rbind(
#'   c(0.682, 1), 
#'   c(0.706, 0.294)
#' )
#' 
#' model_set <- rbind(
#'   c(0.01, 0.9, "weibull"),
#'   c(0.01, 0.99, "lognormal")
#' )
#' 
#' equi <- check_optimality (
#'     best_design = best_design, 
#'     model_set = model_set, 
#'     design_info = design_info
#' )
#'                           
#' print(equi)                          
#' 
#' @references 
#' \enumerate{
#'   \item Müller, C. H., & Pázman, A. (1998). Applications of necessary and sufficient conditions for maximin efficient designs. Metrika, 48, 1–19.
#'   \item Huang, M.-N. L., & Lin, C.-S. (2006). Minimax and maximin efficient designs for estimating the location-shift parameter of parallel models with dual responses. Journal of Multivariate Analysis, 97(1), 198–210.
#' }
#' 
#' @export
check_optimality <- function(best_design, model_set, design_info, seed = 42) {
    
    # check design_info
    if (!inherits(design_info, "DesignInfo")) {
        stop("`design_info` must be a `DesignInfo` object.")
    }
    
    
    # check best_design
    check_design_measure(best_design, design_info)
    
    
    # check model_set
    check_model_set(model_set, design_info)
    
    
    # transform design into particle
    j <- nrow(best_design)
    k <- ncol(best_design)
    transform_prop <- get_transform_prop(best_design[j,])
    best_particle <- c(t(best_design[1:j-1, 1:k]))
    best_particle <- c(best_particle, transform_prop)
    
    
    # coding distribution
    model_set_code <- model_set
    dist_vec <- model_set[, ncol(model_set)]
    model_set_code[, ncol(model_set_code)] <- ifelse(dist_vec == "weibull", 1, 
                                                     ifelse(dist_vec == "lognormal", 2, 3))
    model_set_code <- matrix(
        as.numeric(model_set_code),
        nrow = nrow(model_set_code)
    )
    
    
    ## seed
    seed <- round(seed, digits = 0)
    set.seed(seed)
    
    
    ## run
    optimality <- equivalence_theorem(best_particle, design_info, model_set_code, seed)
    optimality$design <- best_design
    
    class(optimality) <- "OptimalityCheck"
    
    return(optimality)
}


#' Update Optimality Check for Optimal Design
#' 
#' Update optimality check for an OptimalALT object
#'
#' @param design An object of class `OptimalALT`
#' @param check An object of class `OptimalityCheck`
#' @return An `OptimalityALT` object with updated optimality check information
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
#'     sigma=0.6, 
#'     x_l = 0.1, 
#'     x_h = 1
#' )
#' 
#' pso_info <- pso_setting(
#'     n_swarm=5, 
#'     max_iter=10, 
#'     early_stopping=5,
#'     tol=0.0001
#' )
#' 
#' example1_lnorm_minimax <- find_optimal_alt(
#'     design_type="minimax", 
#'     distribution="lognormal",
#'     design_info=design_info, 
#'     pso_info=pso_info, 
#'     coef_lower=c(10^-6, 0.7),
#'     coef_upper=c(10^-4, 0.9), 
#'     highest_level = TRUE,
#'     verbose=TRUE, 
#'     n_threads = 1)
#' 
#' summary(example1_lnorm_minimax)
#' 
#' opt_design <- extract_design(example1_lnorm_minimax)
#' 
#' model_set <- rbind(
#'     c(10^-6, 0.99, "lognormal"),
#'     c(10^-4, 0.99, "lognormal")
#' )
#' 
#' equi <- check_optimality (
#'     best_design=opt_design,
#'     model_set=model_set, 
#'     design_info=design_info
#' )
#' 
#' opt_check_update_design <- update_optimality_check(example1_lnorm_minimax, equi)
#' 
#' summary(opt_check_update_design)
#' 
#' @export 
update_optimality_check <- function(design, check) {
    stopifnot(inherits(design, "OptimalALT"))
    stopifnot(inherits(check, "OptimalityCheck"))
    
    design$max_directional_derivative <- check$max_directional_derivative
    design$model_set <- check$model_set
    design$model_weight <- check$model_weight
    design$equivalence_data <- check$equivalence_data
    
    return(design)
}

