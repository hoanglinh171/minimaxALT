#' @export
set_design_info <- function(k_levels, j_factor, n_unit, censor_time, 
                            p, use_cond, sigma, 
                            x_l = 0, x_h = 1,
                            opt_type = "C", reparam = TRUE,
                            degenerate = FALSE) {
  
  design_info_list <- list()

  design_info_list$n_support = k_levels
  design_info_list$n_factor = j_factor
  design_info_list$n_unit = n_unit
  design_info_list$censor_time = censor_time
  design_info_list$sigma = sigma
  design_info_list$p = p
  design_info_list$use_cond = use_cond
  design_info_list$opt_type = opt_type
  design_info_list$reparam = reparam
  design_info_list$x_l = x_l
  design_info_list$x_h = x_h
  design_info_list$degenerate = degenerate
  
  return(design_info_list)
}


#' @export
pso_setting <- function(n_swarm = 32, max_iter = 100,
                        early_stopping = 10, tol = 0.01, c1 = 2.05, c2 = 2.05,
                        w0 = 1.2, w1 = 0.2, w_var = 0.8, vk = 4) {
  list(
    n_swarm = n_swarm, max_iter = max_iter, 
    early_stopping = early_stopping, tol = tol, c1 = c1, c2 = c2, w0 = w0, w1 = w1, w_var = w_var,
    vk = vk
  )
}


#' @export
initialize_values <- function(init_swarm = NULL,
                              init_local = NULL, init_coef_mat = NULL) {
  list(
    init_swarm = init_swarm, 
    init_local = init_local, init_coef_mat = init_coef_mat
  )
}