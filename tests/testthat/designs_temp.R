devtools::load_all()
library(minimaxALT)

design_info <- set_design_info(k_levels=4, j_factor=1, n_unit=300,
                               censor_time=183, p=0.1, use_cond=c(0),
                               sigma=0.6, x_l = 0.1, x_h = 1)

pso_info <- pso_setting(n_swarm=90, max_iter=180, early_stopping=5, tol=0.0001)

pl_arr <- c(1e-2)
ph_arr <- c(0.9, 0.99, 1)

for (pl in pl_arr) {
  for (ph in ph_arr) {
    res <- find_optimal_alt(design_type=2, distribution=3, design_info=design_info,
                            pso_info=pso_info,
                            # coef=c(1.82*10^-6, 3.17*10^-6, 1),
                            coef_lower = c(1e-5, 0.7),
                            coef_upper = c(pl, ph),
                            # highest_level = FALSE,
                            # init_values = init_values,
                            verbose = TRUE, n_threads = 10
    )
  
    plot(res)
    saveRDS(res, paste0("robust_105", pl, "_07", ph, "_xl01_p01.RDS"))
  }
}
    
    
res <- readRDS("~/minimax_alt_results/minimaxALT/robust_1050.001_071_xl01_p01.RDS")
plot(res)
