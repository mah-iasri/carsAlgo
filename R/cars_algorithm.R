#' Create an object of the CARS algorithm
#'
#' @description
#' The `CARSAlgorithm()` function creates a configuration object for the
#' Competitive Adaptive Reweighted Sampling (CARS) algorithm. Pass this object
#' to \code{\link{fit.CARSAlgorithm}} to run variable selection on your high dimensional dataset.
#'
#' @param max_iter     Maximum number of CARS iterations. Default `100`.
#' @param N            Number of Monte Carlo sub-sampling runs per iteration. Default `50`.
#' @param cv_folds     Number of folds for k-fold cross-validation. Default `5`.
#' @param random_state Integer seed for reproducibility. Default `42`.
#'
#' @return An object of class `"CARSAlgorithm"` - a named list of
#'   hyperparameters to be passed to \code{\link{fit.CARSAlgorithm}}.
#'
#' @examples
#' cars_obj <- CARSAlgorithm(max_iter = 20, N = 30, cv_folds = 5)
#' cars_obj
#'
#' @seealso \code{\link{fit.CARSAlgorithm}}
#' @export
CARSAlgorithm <- function(max_iter     = 100,
                          N            = 50,
                          cv_folds     = 5,
                          random_state = 42) {

  stopifnot(
    is.numeric(max_iter)     && max_iter     >= 1,
    is.numeric(N)            && N            >= 1,
    is.numeric(cv_folds)     && cv_folds     >= 2,
    is.numeric(random_state)
  )

  structure(
    list(
      max_iter     = as.integer(max_iter),
      N            = as.integer(N),
      cv_folds     = as.integer(cv_folds),
      random_state = as.integer(random_state)
    ),
    class = "CARSAlgorithm"
  )
}


#' Print method for CARSAlgorithm objects
#' @param x A `CARSAlgorithm` object.
#' @param ... Ignored.
#' @export
print.CARSAlgorithm <- function(x, ...) {
  cat("CARS Algorithm Configuration\n")
  cat("  max_iter    :", x$max_iter,     "\n")
  cat("  N (MC runs) :", x$N,            "\n")
  cat("  cv_folds    :", x$cv_folds,     "\n")
  cat("  random_state:", x$random_state, "\n")
  invisible(x)
}
