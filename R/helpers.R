
#### Transformation functions --------------------------------------------------

## Transform coefficients by sigmoid function

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

## Get proportion of designs

get_proportion <- function(transform_prop) {
  stopifnot(is.vector(transform_prop),
            all(is.finite(transform_prop)),
            all(transform_prop >= 0),
            all(transform_prop <= 1)
  )
  
  return(transform_proportion(transform_prop))
}

## Get transformed proportion for PSO

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

