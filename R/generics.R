# R/generics.R

#' Fit a Model Object to Data
#'
#' @description
#' Generic function for fitting model objects to data. Methods are
#' dispatched based on the class of `x`.
#'
#' @param cars_obj   A model configuration object (e.g., a `CARSAlgorithm` object).
#' @param ...        Additional arguments passed to the specific method.
#'
#' @return Depends on the method. See \code{\link{fit.CARSAlgorithm}}.
#'
#' @seealso \code{\link{fit.CARSAlgorithm}}
#' @export
fit <- function(cars_obj, ...) {
  UseMethod("fit")
}
