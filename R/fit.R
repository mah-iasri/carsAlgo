#' Fits a CARS Object to any high dimensional dataset
#'
#' @description
#' Applies the CARS algorithm to a high-dimensional data matrix `X` and
#' response vector `y`, iteratively selecting the optimal variable subset
#' via Monte Carlo enabled PLS regression and adaptive reweighted sampling techniques.
#'
#' @param cars_obj       A `CARSAlgorithm` object created by \code{\link{CARSAlgorithm}}.
#' @param X              Numeric matrix of predictors (n_samples x n_features).
#' @param y              Numeric response vector of length n_samples.
#' @param max_components Integer cap on PLS latent components. Default `10`.
#' @param plot           Logical. Whether to display and save the RMSECV curve. Default `TRUE`.
#' @param plot_path      File path for saving the RMSECV plot. Default `"../cars_rmsecv_curve.jpg"`.
#' @param ...            Currently unused.
#'
#' @details
#' This function iteratively:
#' \enumerate{
#'   \item Sub-samples the calibration set (Monte Carlo, `N` runs per iteration).
#'   \item Fits a PLS model and extracts regression coefficients.
#'   \item Selects variables by Adaptive Reweighted Sampling (ARS) proportional
#'         to absolute coefficient magnitude.
#'   \item Evaluates the subset via k-fold cross-validation (RMSECV).
#'   \item Retains the best subset and repeats with an exponentially shrinking
#'         variable set.
#' }
#'
#' @return A named list with:
#' \describe{
#'   \item{`best_features`}{Sorted 1-based column indices of selected features.}
#'   \item{`best_rmsecv`}{Lowest RMSECV achieved across all iterations.}
#'   \item{`rmsecv_history`}{Numeric vector of best RMSECV per iteration.}
#'   \item{`num_features_history`}{Integer vector of feature count per iteration.}
#'   \item{`plot`}{A `ggplot2` object of the RMSECV curve.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(100 * 200), nrow = 100)
#' y <- X[, 5] * 2 + X[, 50] * -1.5 + rnorm(100, sd = 0.5)
#'
#' cars_obj <- CARSAlgorithm(max_iter = 15, N = 30, cv_folds = 5)
#' result   <- fit(cars_obj, X, y, max_components = 8)
#'
#' cat("Best RMSECV      :", result$best_rmsecv, "\n")
#' cat("Selected features:", result$best_features, "\n")
#' }
#'
#' @seealso \code{\link{CARSAlgorithm}}
#'
#' @importFrom pls plsr RMSEP
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal ggsave
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats coef predict
#' @importFrom rlang .data
#' @export
fit.CARSAlgorithm <- function(cars_obj,
                              X,
                              y,
                              max_components = 10L,
                              plot           = TRUE,
                              plot_path      = "../cars_rmsecv_curve.jpg",
                              ...) {

  # Unpack hyperparameters from the CARSAlgorithm object
  max_iter     <- cars_obj$max_iter
  N            <- cars_obj$N
  cv_folds     <- cars_obj$cv_folds
  random_state <- cars_obj$random_state

  # ---- validate inputs ----
  if (!is.matrix(X))  X <- as.matrix(X)
  if (!is.numeric(X)) stop("`X` must be a numeric matrix.")
  if (!is.numeric(y)) stop("`y` must be a numeric vector.")
  if (nrow(X) != length(y))
    stop("`X` and `y` must have the same number of observations.")
  max_components <- as.integer(max_components)

  set.seed(random_state)

  n_features           <- ncol(X)
  best_rmsecv          <- Inf
  best_features        <- NULL
  rmsecv_history       <- numeric(0)
  num_features_history <- integer(0)
  current_features     <- seq_len(n_features)

  message("[cars] Starting CARS - ", n_features, " initial features, ", max_iter, " max iterations, ", N, " MC runs each.")

  pb <- utils::txtProgressBar(min = 0, max = max_iter, style = 3)

  for (iteration in seq_len(max_iter)) {

    n_select <- .exponential_decreasing_function(
      k      = iteration - 1L,
      n_feat = n_features
    )

    if (n_select < 2L) {
      message(sprintf("\n[cars] Stopping at iteration %d: budget < 2 features.", iteration))
      break
    }

    mc_result <- .run_monte_carlo(
      X                = X,
      y                = y,
      current_features = current_features,
      n_select         = n_select,
      max_components   = max_components,
      iteration        = iteration,
      N                = N,
      cv_folds         = cv_folds,
      random_state     = random_state
    )

    if (is.null(mc_result)) {
      message(sprintf("\n[cars] No successful MC runs at iteration %d. Stopping.", iteration))
      break
    }

    rmsecv_history       <- c(rmsecv_history, mc_result$rmsecv)
    num_features_history <- c(num_features_history, length(mc_result$features))

    if (mc_result$rmsecv < best_rmsecv) {
      best_rmsecv   <- mc_result$rmsecv
      best_features <- mc_result$features
    }

    current_features <- mc_result$features
    utils::setTxtProgressBar(pb, iteration)
  }

  close(pb)
  if (!is.null(best_features)) best_features <- sort(best_features)

  #### Generate plots
  best_idx        <- which.min(rmsecv_history)
  best_feature_pt <- num_features_history[best_idx]

  plot_df <- data.frame(
    num_features = num_features_history,
    rmsecv       = rmsecv_history
  )

  # Best point as separate data frame — columns named to avoid
  # clashing with function parameters x and y
  best_pt_df <- data.frame(
    pt_features = best_feature_pt,
    pt_rmsecv   = best_rmsecv
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$num_features, y = .data$rmsecv)   # .data$ fixes num_features and rmsecv
  ) +
    ggplot2::geom_line(color = "black", linewidth = 0.8) +
    ggplot2::geom_point(color = "black", size = 2) +
    ggplot2::geom_point(
      data = best_pt_df,
      ggplot2::aes(x = .data$pt_features, y = .data$pt_rmsecv), # .data$ fixes x
      color      = "darkred",
      shape      = 15,
      size       = 3,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x     = "Number of variables",
      y     = "RMSECV",
      title = paste0("Selected variables: ", length(best_features),
                     "  (RMSECV: ", round(best_rmsecv, 4), ")")
    ) +
    ggplot2::theme_minimal()

  print(p)

  if (plot) {
    dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = plot_path, plot = p, width = 10, height = 6, dpi = 300)
    message("Plot saved to: ", plot_path)
  }

  list(
    best_features        = best_features,
    best_rmsecv          = best_rmsecv,
    rmsecv_history       = rmsecv_history,
    num_features_history = num_features_history,
    plot                 = p
  )
}
