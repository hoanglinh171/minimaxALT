#' Extract design from an OptimalALT object
#'
#' @param object An object of class `OptimalALT`.
#' @return The extracted design.
#' @export
extract <- function(object) {
    UseMethod("extract")
}


#' Update optimality check for an OptimalALT object
#'
#' @param object An object of class `OptimalALT`
#' @param check An object of class `OptimalityCheck`
#' @return An `OptimalityALT` object with updated optimality check information
#' @export
update_optimality_check <- function(object, check) {
    UseMethod("update_optimality_check")
}
