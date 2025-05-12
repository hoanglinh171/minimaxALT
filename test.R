library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(ggplot2)
library(parallel)

library(devtools)

use_git()

Rcpp::compileAttributes()
# devtools::clean_dll()
devtools::load_all()

# Sys.setenv(PKG_CPPFLAGS = "-g -O0")  # -g: debug info, -O0: no optimization
# devtools::install("C:/libs/minimaxALT")
# devtools::install("C:/libs/minimaxALT",
#                   force = TRUE,
#                   args = "--with-keep.source --debug")
# 
# devtools::uninstall("C:/libs/minimaxALT")


library(minimaxALT)


coef = c(0.001, 0.9)

coef_lower = c(1e-6, 0.7)
coef_upper = c(1e-2, 1)
n_support = 4
n_factor = 1
n_unit = 300
censor_time = 183
sigma = 0.6
p = 0.1
use_cond = c(0)
distribution = 1  # 1: weibull, 2: log normal, 3: model robust
design_type = 1 # 1: locally optimal, 2: minimax


get_outbound_sigmoid(c(1e-6, 5e-3, 1), c(1e-6, 1e-6, 0.7), c(1e-2, 1e-2, 1))

init_coef <- rbind(
  c(1e-2, 1e-2, 1),
  c(1e-6, 5e-3, 1),
  c(1e-2, 1e-2, 0.7),
  c(1e-6, 5e-3, 0.7),
  c(1e-6, 1e-6, 1),
  c(1e-6, 1e-6, 0.7)
)

init_values = initialize_values(init_coef_mat = init_coef)

design_info <- set_design_info(k_levels=4, j_factor=2, n_unit=170, 
                               censor_time=1000, p=0.001, use_cond=c(0, 0), 
                               sigma=0.6743, x_l = 0, x_h = 1)

pso_info <- pso_setting(n_swarm=256, max_iter=512, early_stopping=10, tol=0.01)

# set.seed(30)
# res <- find_optimal_alt(design_type=2, distribution=3, design_info=design_info, 
#                         pso_info=pso_info, 
#                         # coef=c(1.82*10^-6, 3.17*10^-6, 1), 
#                         coef_lower = c(1e-6, 0.7),
#                         coef_upper = c(1e-2, 0.99),
#                         highest_level = FALSE,
#                         verbose = TRUE, n_threads = 10)

set.seed(30)

# library(pryr)
mem_change(res <- find_optimal_alt(design_type=2, distribution=1, design_info=design_info, 
                        pso_info=pso_info,
                        # coef=c(1.82*10^-6, 3.17*10^-6, 1),
                        coef_lower = c(1e-6, 1e-6, 0.7),
                        coef_upper = c(1e-2, 1e-2, 1),
                        # highest_level = FALSE,
                        init_values = init_values,
                        verbose = TRUE, n_threads = 16
                        )
)


plot(res, x_l = 0, x_h =1, nlevels = 30)


saveRDS(res, "1f_robust_1060.0001_0.10.7_xl0_p0.1.txt")

summary(res)
 ### ---------------------------------------------------------------------------
### Equivalence theorem

model_set = rbind(
  c(1e-2, 0.9, 2),
  c(1e-6, 0.9, 1),
  c(1e-6, 0.9, 2),
  c(1e-2, 0.9, 1),
  c(1e-6, 0.7, 1),
  c(1e-6, 0.7, 2)
)


equi_check <- check_equivalence_theorem (best_particle = as.vector(res$g_best), model_set = model_set,
                                      design_info = design_info)


equi_check$max_directional_derivative
equi_check$model_weight

equi <- as.data.frame(equi_check$equivalence_data)

colnames(equi) <- c("Stress level", "Directional derivative", "Point")
points <- equi[equi$Point == 1, ]
# points <- points[c(1, 2, 4), ]

ggplot(equi, aes(x=`Stress level`, y=`Directional derivative`)) +
  geom_hline(yintercept = 1, color="darkgrey") +
  geom_line() +
  annotate("point", x = points$`Stress level`, y = points$`Directional derivative`, colour = "blue") +
  xlim(0, 1) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"))


# equi_106105_0709_xl0


fisher <- rbind(
  c(681.85381879022, 158.07433237550, 96.48949321427, 156.32195374133),
  c(158.07433237550, 373.88891312666, 228.22390851575, 369.74405979651),
  c(96.48949321427, 228.22390851575, 185.78348706880, 225.59422280774),
  c(156.32195374133, 369.74405979651, 225.59422280774, 365.64506055108)
)

rank(fisher)

solve(fisher)
