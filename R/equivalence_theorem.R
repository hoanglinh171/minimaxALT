#' Check Equivalence Theorem for Optimal Design
#'
#' Evaluates whether a design satisfies the equivalence theorem.
#'
#' @param best_particle A vector containing the best particle's position (i.e., stress levels and transformed allocated proportion).
#' @param model_set A matrix of models, including parameters and distribution, that maximize the optimality criteria with the given best particle's position.
#' @param design_info A list containing design parameters such as factor levels, number of units, and other settings.
#'
#' @return A list containing
#' \describe{
#' \item{max_directional_derivative}{Maximum directional derivative within design space.}
#' \item{model_set}{The model set that is input.}
#' \item{model_weight}{The weight assigned to each model in the model set.}
#' \item{equivalence_data}{Generated designs and their corresponding directional derivative given the optimal design \code{best_particle}. Each design is a combination of factors with value in [0, 1]. These designs are data for plotting equivalence theorem plot.}
#' }
#' @examples
#' three-dimensional function
#' objf <- function(x, loc) {
#'   val <- 0
#'   for (i in 1:length(x)) {
#'     val <- val + (x[i] - loc)^2
#'   }
#'   return(val)
#' }
#' 
#' upp_bound <- rep(5, 3)
#' low_bound <- rep(-5, 3)
#' loc_shift <- 1
#'
#' alg_setting <- getPSOInfo(nSwarm = 32, maxIter = 100, psoType = "basic")
#' res <- globpso(objFunc = objf, lower = low_bound, upper = upp_bound, 
#'                PSO_INFO = alg_setting, loc = loc_shift)
#' res$par
#' res$val
#' @references 
#' \enumerate{
#'   \item Bonyadi, M. R., & Michalewicz, Z. (2014). A locally convergent rotationally invariant particle swarm optimization algorithm. Swarm intelligence, 8(3), 159-198.
#'   \item Cheng, R., & Jin, Y. (2014). A competitive swarm optimizer for large scale optimization. IEEE transactions on cybernetics, 45(2), 191-204.
#   \item Eberhart, R. & Kennedy, J. (1995). A new optimizer using particle swarm theory. In The 6th International Symposium on Micro Machine and Human Science, pages 39-43. IEEE.
#'   \item Shi, Y., & Eberhart, R. (1998, May). A modified particle swarm optimizer. In Evolutionary Computation Proceedings, 1998. IEEE World Congress on Computational Intelligence., The 1998 IEEE International Conference on (pp. 69-73). IEEE.
#'   \item Stehlík, M., Chen, P. Y., Wong, W. K., and Kiseľák, J. (2024). A double exponential particle swarm optimization with non-uniform variates as stochastic tuning and guaranteed convergence to a global optimum with sample applications to finding optimal exact designs in biostatistics. Applied Soft Computing, 163, 111913.
#'   \item Sun, J., Feng, B., and Xu, W. (2004a). Particle swarm optimization with particles having quantum behavior. In Evolutionary Computation, 2004. CEC2004. Congress on, volume 1, pages 325-331. IEEE.
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

