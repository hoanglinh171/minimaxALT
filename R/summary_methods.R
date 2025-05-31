
#' @export
summary.OptimalALT <- function(object, ...) {
  stopifnot(inherits(object, "OptimalALT"))
  
  n_factor <- length(object$coef_best) - 1
  n_support <- (length(object$g_best) + 1) / (n_factor + 1)
  
  stress_levels <- matrix(object$g_best[1:(n_support*n_factor)], 
                          ncol = n_support, byrow=TRUE)
  
  prop <- get_proportion(object$g_best[(n_support*n_factor + 1):length(object$g_best)])
  
  design <- rbind(stress_levels, prop)
  
  level_names <- paste0("X", 1:n_factor)
  level_names <- c(level_names, "W")
  rownames(design) <- level_names
  
  cat("Summary of generated optimal ALT design\n")
  cat("-----------------------------------------------\n")
  cat("X: Stress levels, W: Corresponding proportion\n")
  print(design)
  cat("\nObjective Value:", object$fg_best, "\n")
  cat("Max directional derivative:", object$max_directional_derivative)

  invisible(object)
}


#' @export
plot.OptimalALT <- function(x, ...) {
  stopifnot(inherits(x, "OptimalALT"))
  
  n_factor <- length(x$coef_best) - 1
  n_support <- (length(x$g_best) + 1) / (n_factor + 1)
  
  stress_levels <- matrix(x$g_best[1:(n_support*n_factor)], 
                          ncol = n_support, byrow=TRUE)
  
  prop <- get_proportion(x$g_best[(n_support*n_factor + 1):length(x$g_best)])
  
  args <- list(...)
  x_l <- ifelse(is.null(args$x_l), 0, args$x_l)
  x_h <- ifelse(is.null(args$x_h), 1, args$x_h)
  nlevels <- ifelse(is.null(args$nlevels), 10, args$nlevels)

  
  if (n_factor == 1) {
    plot_one_factor(x$equivalence_data, proportion = prop, x_l=x_l, x_h=x_h)
    
  } else if (n_factor == 2) {
    plot_two_factor(x$equivalence_data, proportion = prop, x_l=x_l, x_h=x_h,
                    nlevels = nlevels)
    
  } else {
    stop("Do not support plotting for ALT with more than 2 factors.")
  }
  
  invisible(x)
}


plot_one_factor <- function(equivalence_data, proportion, x_l, x_h) {
  equi <- as.data.frame(equivalence_data)
  
  colnames(equi) <- c("Stress level", "Directional derivative", "Point")
  points <- equi[equi$Point == 1, ]
  
  valid_idx <- proportion >= 0.001
  points <- points[valid_idx, ]
  
  
  p <- ggplot(equi, aes(x=`Stress level`, y=`Directional derivative`)) +
    geom_hline(yintercept = 1, color="darkgrey") +
    # geom_vline(xintercept = x_l, color="red", linetype="dashed") +
    geom_vline(xintercept = x_h, color="red", linetype="dashed") +
    geom_line() +
    annotate("point", x = points$`Stress level`, y = points$`Directional derivative`, colour = "blue") +
    xlim(0, 1) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.line = element_line(color = "black"))
  
  print(p)
  
}


plot_two_factor <- function(equivalence_data, proportion, x_l, x_h, nlevels) {
  equi <- as.data.frame(equivalence_data)
  colnames(equi) <- c("x1", "x2", "dd", "Point")
  points <- equi[equi$Point == 1, ]
  valid_idx <- proportion >= 0.001
  points <- points[valid_idx, ]
  points$dd <- round(points$dd, digits = 2)
  equi <- equi[equi$Point == 0, ]
  
  
  x_vals <- sort(unique(equi$x1))
  y_vals <- sort(unique(equi$x2))
  
  z <- matrix(equi$dd, ncol = length(y_vals), nrow = length(x_vals), byrow = TRUE)
  
  contour(x_vals, y_vals, z, nlevels = nlevels,
          xlab = "Stress level 1", ylab = "Stress level 2")
  
  rect(xleft = x_l, ybottom = x_l, xright = x_h, ytop = x_h, border = "red", lty = "dashed")
  points(points$x1, points$x2, col = "blue", pch = 19, cex = 0.7)
  text(points$x1, points$x2 + 0.03, points$dd,col = "blue", cex = 0.7)
  
}
