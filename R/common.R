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

get_proportion <- function(dirichlet_prop) {
  stopifnot(is.vector(dirichlet_prop),
            all(is.finite(dirichlet_prop)),
            all(dirichlet_prop >= 0),
            all(dirichlet_prop <= 1)
  )
  
  return(transform_proportion(dirichlet_prop))
}
