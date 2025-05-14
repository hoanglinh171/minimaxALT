#' Check Equivalence Theorem for Optimal Design
#'
#' Evaluates whether a design satisfies the equivalence theorem.
#'
#' @param best_particle A vector containing the best particle's position (i.e., stress levels and transformed allocated proportion).
#' @param model_set A matrix of models, including parameters and distribution, that maximize the optimality criteria with the given best particle's position.
#' @param design_info A list containing design parameters such as factor levels, number of units, and other settings.
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
#' best_particle <- c(0.682, 1, 0.706)
#' 
#' model_set <- rbind(
#'   c(0.01, 0.9, 1),
#'   c(0.01, 0.99, 2))
#' 
#' equi <- check_equivalence_theorem (best_particle=best_particle, 
#'                                     model_set=model_set, 
#'                                     design_info=design_info)
#' 
#' equi$max_directional_derivative
#' 
#' @references 
#' \enumerate{
#'   \item Müller, C. H., & Pázman, A. (1998). Applications of necessary and sufficient conditions for maximin efficient designs. Metrika, 48, 1–19.
#'   \item Huang, M.-N. L., & Lin, C.-S. (2006). Minimax and maximin efficient designs for estimating the location-shift parameter of parallel models with dual responses. Journal of Multivariate Analysis, 97(1), 198–210.
#' }
#' @name check_equivalence_theorem
#' @rdname check_equivalence_theorem
#' @importFrom Rcpp evalCpp cppFunction sourceCpp
#' @export
check_equivalence_theorem <- function(best_particle, model_set, design_info) {
  
  stopifnot(design_info$opt_type == "C" || design_info$opt_type == "D")
  
  stopifnot(is.numeric(design_info$n_support), is.numeric(design_info$n_factor), 
            is.numeric(design_info$n_unit), 
            is.numeric(design_info$censor_time), is.numeric(design_info$sigma), 
            is.numeric(design_info$p), 
            is.numeric(design_info$x_l), is.numeric(design_info$x_h))
  
  stopifnot(is.logical(design_info$reparam), is.logical(design_info$degenerate))
  
  use_cond = c(design_info$use_cond)
  stopifnot(design_info$n_factor == length(use_cond))
  
  stopifnot(is.vector(best_particle))
  stopifnot(is.matrix(model_set))
  
  ## Define design info
  design_info$use_cond = use_cond
  
  return(equivalence_theorem(best_particle, design_info, model_set))
}

