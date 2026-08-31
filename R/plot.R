#' Plot an OptimalALT Object
#'
#' Plot the verification of design optimality.
#'
#' @param x An object of class `OptimalALT`, typically returned by
#'   \code{minimax_alt()}.
#' @param x_l Numeric. Lower bound of the stress range. Default is \code{0}.
#' @param x_h Numeric. Upper bound of the stress range. Default is \code{1}.
#' @param nlevels Integer. Number of grid levels used for plotting the 
#' optimality check for a two-factor design. Default is
#'   \code{10}.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns the original `OptimalALT` object.
#'
#' @examples
#' \dontrun{
#' # Suppose `result` is an OptimalALT object returned by find_optimal_alt().
#' plot(result)
#'
#' # Specify the stress range explicitly.
#' plot(result, x_l = 0, x_h = 1)
#'
#' # For a two-factor design, increase the plotting grid resolution.
#' plot(result, x_l = 0, x_h = 1, nlevels = 20)
#' }
#'
#' @export
plot.OptimalALT <- function(x, x_l = 0, x_h = 1, nlevels = 10, ...) {
    stopifnot(inherits(x, "OptimalALT"))
    
    n_factor <- length(x$coef_best) - 1
    n_support <- (length(x$g_best) + 1) / (n_factor + 1)
    
    prop <- get_proportion(x$g_best[(n_support*n_factor + 1):length(x$g_best)])
    
    x_l <- ifelse(is.null(x_l), 0, x_l)
    x_h <- ifelse(is.null(x_h), 1, x_h)
    nlevels <- ifelse(is.null(nlevels), 10, nlevels)
    
    
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


## plot.OptimalALT() helper functions ------------------------------------------

utils::globalVariables(c("stress_level", "dir_deriv"))

plot_one_factor <- function(equivalence_data, proportion, x_l, x_h) {
    equi <- as.data.frame(equivalence_data)
    
    colnames(equi) <- c("stress_level", "dir_deriv", "point")
    points <- equi[equi$point == 1, ]
    
    valid_idx <- proportion >= 0.001
    points <- points[valid_idx, ]
    
    
    p <- ggplot(equi, aes(x=stress_level, y=dir_deriv)) +
        geom_hline(yintercept = 1, color="darkgrey") +
        geom_vline(xintercept = x_l, color="red", linetype="dashed") +
        geom_vline(xintercept = x_h, color="red", linetype="dashed") +
        geom_line() +
        annotate("point", x = points$stress_level, y = points$dir_deriv, colour = "blue") +
        xlim(0, 1) +
        xlab("Stress Level") +
        ylab("Directional Derivative") +
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
