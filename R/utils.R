# ============================================================================ #
#  utils.R - Internal helper functions for the CARS algorithm                 #
#  All functions are prefixed with `.` to signal they are private.            #
#  None are exported; accessible internally via the package namespace.        #
# ============================================================================ #


# ---------------------------------------------------------------------------- #
#  .exponential_decreasing_function                                            #
# ---------------------------------------------------------------------------- #

#' Exponential Decreasing Schedule (Internal)
#'
#' Computes the number of variables to retain at CARS iteration `k` using
#' the schedule: \eqn{\lceil (n\_feat / 2^b) \cdot e^{-k/b} \rceil}.
#'
#' @param k      Integer. Current iteration index (0-based).
#' @param n_feat Integer. Total number of features at the start of CARS.
#' @param b      Numeric. Decay rate parameter. Default `2`.
#'
#' @return Integer. Number of features to retain at this iteration.
#' @keywords internal
.exponential_decreasing_function <- function(k, n_feat, b = 2) {
  a <- n_feat / (2 ^ b)
  as.integer(ceiling(a * exp(-k / b)))
}


# ---------------------------------------------------------------------------- #
#  .adaptive_reweighted_sampling                                               #
# ---------------------------------------------------------------------------- #

#' Adaptive Reweighted Sampling (Internal)
#'
#' Samples `num_features` variable indices without replacement, with
#' probabilities proportional to the absolute values of PLS regression
#' coefficients. Falls back to uniform probabilities if weights are all
#' zero, non-finite, or contain `NA`.
#'
#' @param weights      Numeric vector of PLS regression coefficients.
#' @param num_features Integer. Number of variable indices to select.
#'
#' @return Integer vector of selected indices (positions within `weights`).
#' @keywords internal
.adaptive_reweighted_sampling <- function(weights, num_features) {
  weights_abs <- abs(weights)

  if (sum(weights_abs) == 0 ||
      any(is.na(weights_abs)) ||
      any(is.infinite(weights_abs))) {
    weights_norm <- rep(1 / length(weights_abs), length(weights_abs))
  } else {
    weights_norm <- weights_abs / sum(weights_abs)
  }

  sample(
    seq_along(weights),
    size    = num_features,
    replace = FALSE,
    prob    = weights_norm
  )
}


# ---------------------------------------------------------------------------- #
#  .get_optimal_components                                                     #
# ---------------------------------------------------------------------------- #

#' Determine Safe Number of PLS Components (Internal)
#'
#' Returns the largest safe number of PLS latent components given the
#' current sample size and feature count, capped at `max_components`.
#' Ensures the value is at least 1 to avoid degenerate PLS models.
#'
#' @param n_samples      Integer. Number of available (calibration) samples.
#' @param n_features     Integer. Number of features in the current subset.
#' @param max_components Integer. Hard upper bound supplied by the user.
#'
#' @return Integer. Number of PLS components to use.
#' @keywords internal
.get_optimal_components <- function(n_samples, n_features, max_components) {
  max_allowed <- min(n_samples - 1L, n_features, max_components)
  max(1L, as.integer(max_allowed))
}


# ---------------------------------------------------------------------------- #
#  .pls_rmsecv                                                                 #
# ---------------------------------------------------------------------------- #

#' Cross-Validated RMSECV via Manual K-Fold (Internal)
#'
#' Fits a PLS model using manual k-fold cross-validation and returns the
#' Root Mean Square Error of Cross-Validation (RMSECV). Uses
#' \code{\link[pls]{plsr}} with the \code{"kernelpls"} algorithm.
#' Returns `Inf` if all folds fail or produce no predictions.
#'
#' @param X_sel   Numeric matrix. Feature matrix for all samples
#'   (n_samples x selected features).
#' @param y       Numeric vector. Response variable (length n_samples).
#' @param ncomp   Integer. Number of PLS latent components to use.
#' @param n_folds Integer. Number of cross-validation folds.
#' @param seed    Integer or `NULL`. Optional seed for reproducible fold
#'   assignment. Default `NULL`.
#'
#' @return Numeric scalar. RMSECV value; `Inf` if computation fails.
#'
#' @importFrom pls plsr
#' @importFrom stats predict
#' @keywords internal
.pls_rmsecv <- function(X_sel, y, ncomp, n_folds, seed = NULL) {

  n <- nrow(X_sel)

  if (!is.null(seed)) set.seed(seed)
  fold_ids  <- sample(rep(seq_len(n_folds), length.out = n))
  sq_errors <- numeric(0)

  for (fold in seq_len(n_folds)) {

    test_idx  <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)

    if (length(train_idx) < 3L || length(test_idx) < 1L) next

    X_train <- X_sel[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test  <- X_sel[test_idx,  , drop = FALSE]
    y_test  <- y[test_idx]

    tryCatch({
      model <- pls::plsr(
        y ~ .,
        data   = data.frame(y = y_train, X_train),
        ncomp  = ncomp,
        scale  = FALSE,
        method = "kernelpls"
      )

      preds <- stats::predict(
        model,
        newdata = data.frame(X_test),
        ncomp   = ncomp
      )[, 1L, 1L]

      sq_errors <- c(sq_errors, (y_test - preds)^2)

    }, error = function(e) NULL)
  }

  if (length(sq_errors) == 0L) return(Inf)
  sqrt(mean(sq_errors))
}


# ---------------------------------------------------------------------------- #
#  .run_monte_carlo                                                            #
# ---------------------------------------------------------------------------- #

#' Monte Carlo Runs for One CARS Iteration (Internal)
#'
#' Executes all `N` Monte Carlo sub-sampling runs for a single CARS iteration.
#' In each run:
#' \enumerate{
#'   \item A random 80\% calibration subset is drawn.
#'   \item A PLS model is fitted and regression coefficients are extracted.
#'   \item Feature indices are selected via Adaptive Reweighted Sampling (ARS).
#'   \item The selected subset is evaluated via k-fold RMSECV on the full data.
#' }
#' Returns the feature subset and RMSECV of the best-performing run, or
#' `NULL` if every run fails.
#'
#' @param X                Numeric matrix. Full predictor matrix (n_samples x n_features).
#' @param y                Numeric vector. Response variable (length n_samples).
#' @param current_features Integer vector. Active feature indices (1-based) for
#'   this iteration.
#' @param n_select         Integer. Target number of features to select.
#' @param max_components   Integer. PLS component cap (from user).
#' @param iteration        Integer. Current CARS iteration index (used for
#'   deterministic per-run seed offsets).
#' @param N                Integer. Number of Monte Carlo runs to execute.
#' @param cv_folds         Integer. Number of cross-validation folds.
#' @param random_state     Integer. Base random seed from the `CARSAlgorithm` object.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`features`}{Integer vector of selected feature indices from the best run.}
#'   \item{`rmsecv`}{Numeric. RMSECV of the best run.}
#' }
#' Returns `NULL` if no run succeeds.
#'
#' @importFrom pls plsr
#' @importFrom stats coef
#' @keywords internal
.run_monte_carlo <- function(X, y,
                             current_features,
                             n_select,
                             max_components,
                             iteration,
                             N,
                             cv_folds,
                             random_state) {

  n_samples          <- nrow(X)
  n_current_features <- length(current_features)
  iteration_rmsecv   <- numeric(0)
  iteration_features <- list()

  for (run in seq_len(N)) {

    # Sub-sample 80% of data as calibration set
    cal_size   <- min(max(10L, as.integer(0.8 * n_samples)), n_samples)
    sample_idx <- sample(n_samples, size = cal_size, replace = FALSE)

    X_cal <- X[sample_idx, current_features, drop = FALSE]
    y_cal <- y[sample_idx]

    ncomp_cal <- .get_optimal_components(cal_size, n_current_features,
                                         max_components)

    tryCatch({

      # Fit PLS on calibration subset
      pls_model <- pls::plsr(
        y ~ .,
        data   = data.frame(y = y_cal, X_cal),
        ncomp  = ncomp_cal,
        scale  = FALSE,
        method = "kernelpls"
      )

      coef_weights <- as.numeric(
        stats::coef(pls_model, ncomp = ncomp_cal, intercept = FALSE)
      )

      # Feature selection via Adaptive Reweighted Sampling
      selected_features <-
        if (n_current_features > n_select) {
          sel_idx <- .adaptive_reweighted_sampling(coef_weights, n_select)
          current_features[sel_idx]
        } else {
          current_features
        }

      # Evaluate selected subset with k-fold RMSECV on the full dataset
      X_selected <- X[, selected_features, drop = FALSE]
      ncomp_cv   <- .get_optimal_components(
        n_samples, length(selected_features), max_components
      )

      rmsecv_val <- .pls_rmsecv(
        X_sel   = X_selected,
        y       = y,
        ncomp   = ncomp_cv,
        n_folds = min(cv_folds, n_samples %/% 2L),
        seed    = random_state + iteration * 1000L + run
      )

      iteration_rmsecv   <- c(iteration_rmsecv, rmsecv_val)
      iteration_features <- c(iteration_features, list(selected_features))

    }, error = function(e) {
      message(sprintf("  [cars] Warning - Iter %d, Run %d: %s",
                      iteration, run, conditionMessage(e)))
    })
  }

  # Return NULL if every run failed
  if (length(iteration_rmsecv) == 0L) return(NULL)

  best_idx <- which.min(iteration_rmsecv)
  list(
    features = iteration_features[[best_idx]],
    rmsecv   = iteration_rmsecv[[best_idx]]
  )
}
