get_outbound_sigmoid <- function(coef_vec, coef_lower, coef_upper) {
  stopifnot(is.vector(coef_vec),
            all(is.finite(coef_vec)),
            all(is.finite(coef_lower)),
            all(is.finite(coef_upper)),
            (length(coef_vec) == length(coef_lower)),
            all(coef_vec >= coef_lower),
            all(coef_vec <= coef_upper)
            )
  
  return(transform_sigmoid(coef_vec, coef_lower, coef_upper))
}

get_proportion <- function(transform_prop) {
  stopifnot(is.vector(transform_prop),
            all(is.finite(transform_prop)),
            all(transform_prop >= 0),
            all(transform_prop <= 1)
  )
  
  return(transform_proportion(transform_prop))
}

get_transform_prop <- function(prop) {
  stopifnot(is.vector(prop),
            all(is.finite(prop)),
            all(prop >= 0),
            all(prop <= 1),
            sum(prop) == 1
  )
  
  transform_prop = prop[1:(length(prop) - 1)]
  if (length(transform_prop) > 1) {
    denom = 1
    for (i in 2:length(transform_prop)) {
      denom = denom - prop[i-1]
      transform_prop[i] = prop[i] / denom 
    }
  }
  
  return(transform_prop)
}

extract_design <- function(object) {
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
    
    return(design)
    
}
