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

