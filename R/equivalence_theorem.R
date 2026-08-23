#' Check Optimality for Optimal Design
#'
#' Evaluates whether a design is optimal by the equivalence theorem.
#'
#' @param best_design A matrix containing stress levels and allocated proportion of the design.
#' @param model_set A matrix of models, including parameters and distribution, that maximize the optimality criteria with the given best particle's position.
#' @param design_info A list containing design parameters such as factor levels, number of units, and other settings.
#' @param seed Seed for reproducibility
#'
#' @return
#' \describe{
#' \item{max_directional_derivative}{Maximum directional derivative within design space.}
#' \item{model_set}{The model set that is input.}
#' \item{model_weight}{The weight assigned to each model in the model set.}
#' \item{equivalence_data}{Generated designs and their corresponding directional derivative given the optimal design \code{best_particle}. Each design is a combination of factors with value in [0, 1]. These designs are data for plotting equivalence theorem plot.}
#' }
#' @examples
#' design_info <- set_design_info(k_levels=2, j_factor=1, n_unit=300, 
#'                                censor_time=183, p=0.1, use_cond=0, sigma=0.6)
#'                                
#' best_design <- rbind(
#'   c(0.682, 1), 
#'   c(0.706, 0.294)
#' )
#' 
#' model_set <- rbind(
#'   c(0.01, 0.9, 1),
#'   c(0.01, 0.99, 2))
#' 
#' equi <- check_optimality (best_design=best_design, 
#'                           model_set=model_set, 
#'                           design_info=design_info)
#' print(equi)                          
#' 
#' @references 
#' \enumerate{
#'   \item Müller, C. H., & Pázman, A. (1998). Applications of necessary and sufficient conditions for maximin efficient designs. Metrika, 48, 1–19.
#'   \item Huang, M.-N. L., & Lin, C.-S. (2006). Minimax and maximin efficient designs for estimating the location-shift parameter of parallel models with dual responses. Journal of Multivariate Analysis, 97(1), 198–210.
#' }
#' @name check_optimality
#' @rdname check_optimality
#' @importFrom Rcpp evalCpp cppFunction sourceCpp
#' @export
check_optimality <- function(best_design, model_set, design_info, seed = 42) {
  
  # transform design into particle
  stopifnot(is.matrix(best_design))
  j <- nrow(best_design)
  k <- ncol(best_design)
  stopifnot(design_info$n_factor == j - 1)
  stopifnot(design_info$n_support == k)
  stopifnot(all(best_design[j,] >= 0),
            all(best_design[j,] <= 1)
            )
  stopifnot(sum(best_design[j,]) == 1)
  
  transform_prop <- get_transform_prop(best_design[j,])
  best_particle <- c(t(best_design[1:j-1, 1:k]))
  best_particle <- c(best_particle, transform_prop)
  
  stopifnot(is.numeric(design_info$n_support), is.numeric(design_info$n_factor), 
            is.numeric(design_info$n_unit), 
            is.numeric(design_info$censor_time), is.numeric(design_info$sigma), 
            is.numeric(design_info$p), 
            is.numeric(design_info$x_l), is.numeric(design_info$x_h))
  
  stopifnot(is.logical(design_info$reparam))
  
  use_cond = c(design_info$use_cond)
  stopifnot(design_info$n_factor == length(use_cond))
  
  stopifnot(is.matrix(model_set))
  
  seed = round(seed, digits = 0)
  
  ## Define design info
  design_info$use_cond = use_cond
  # stopifnot(design_info$opt_type == "C")
  design_info$opt_type = "C"
  
  optimality <- equivalence_theorem(best_particle, design_info, model_set, seed)
  optimality$design <- best_design
  
  class(optimality) <- "OptimalityCheck"
  
  return(optimality)
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
}

