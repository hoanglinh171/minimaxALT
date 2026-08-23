design_info <- set_design_info(k_levels=3, j_factor=1, n_unit=300,
                               censor_time=183, p=0.1, use_cond=0, sigma=0.6)

best_design <- rbind(
  c(1, 0.488, 0.828),
  c(0.260, 0.394, 0.346)
)

model_set <- rbind(
  c(10^-6, 0.99, 1),
  c(10^-4, 0.99, 2)
)

equi <- check_equivalence_theorem (best_design=best_design,
                                   model_set=model_set, 
                                   design_info=design_info)

equi$max_directional_derivative