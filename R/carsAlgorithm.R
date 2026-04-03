#' Competitive Adaptive Reweighted Sampling (CARS) Algorithm
#'
#' @description
#' Creates a CARS algorithm object that exposes a `$fit()` method.
#' The constructor captures all tuning parameters in a closure; private
#' helper functions are nested inside and never exported.
#'
#' CARS iteratively:
#' \enumerate{
#'   \item Sub-samples the calibration set (Monte Carlo, `N` runs per iteration).
#'   \item Fits a PLS model and extracts regression coefficients.
#'   \item Selects variables by Adaptive Reweighted Sampling (ARS) proportional to absolute coefficient magnitude.
#'   \item Evaluates the subset via k-fold cross-validation (RMSECV).
#'   \item Retains the best subset and repeats with an exponentially shrinking
#'         variable budget.
#' }
#'
#' @param max_iter     Maximum number of CARS iterations. Default `100`.
#' @param N            Number of Monte Carlo sub-sampling runs per iteration.  Default `50`.
#' @param cv_folds     Number of folds for k-fold cross-validation. Default `5`.
#' @param random_state Integer seed for reproducibility. Default `42`.
#'
#' @return A named list with one public element:
#' \describe{
#'   \item{`fit(X, y, max_components = 10,  save_plot= TRUE, plot_path= "./cars_rmsecv_curve.jpg")`}{Fits the CARS object on high dimentsional data matrix `X`
#'     and response vector `y`. Returns a named list with:
#'     \itemize{
#'       \item `best_features` — sorted 1-based column indices of selected features.
#'       \item `best_rmsecv`      — lowest RMSECV achieved across all iterations.
#'       \item `rmsecv_history`   — numeric vector of best RMSECV per iteration.
#'       \item `num_features_history` — integer vector of feature count per iteration.
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate high dimensional data
#' set.seed(1)
#' X <- matrix(rnorm(100 * 200), nrow = 100)
#' y <- X[, 5] * 2 + X[, 50] * -1.5 + rnorm(100, sd = 0.5)
#'
#' # Instantiate and fit
#' cars <- CARSAlgorithm(max_iter = 15, N = 30, cv_folds = 5, random_state = 42)
#' result <- cars$fit(X, y, max_components = 8)
#'
#' cat("Best RMSECV      :", result$best_rmsecv, "\n")
#' cat("Selected features:", result$best_features, "\n")
#' }
#'
#'
# Automatically install required packages if not installed
#' @importFrom pls plsr
#' @importFrom ggplot2 ggplot
#' @export
CARSAlgorithm <- function(max_iter     = 100,
                          N            = 50,
                          cv_folds     = 5,
                          random_state = 42) {

  # ---- input validation ----
  stopifnot(
    is.numeric(max_iter)     && max_iter     >= 1,
    is.numeric(N)            && N            >= 1,
    is.numeric(cv_folds)     && cv_folds     >= 2,
    is.numeric(random_state)
  )
  max_iter     <- as.integer(max_iter)
  N            <- as.integer(N)
  cv_folds     <- as.integer(cv_folds)
  random_state <- as.integer(random_state)


  # ------------------------------------------------------------------ #
  #  PRIVATE helpers (enclosed — not exported, not visible to users)    #
  # ------------------------------------------------------------------ #

  # ------------------------------------------------------------------
  # .exponential_decreasing_function
  # Returns the number of variables to keep at iteration k using an
  # exponential decay schedule:  ceil( (n_feat / 2^b) * exp(-k / b) )
  #
  # @param k       Iteration index, 0-based.
  # @param n_feat  Number of variables at the start of CARS.
  # @param b       Decay rate parameter (default 2).
  # ------------------------------------------------------------------
  .exponential_decreasing_function <- function(k, n_feat, b = 2) {
    a <- n_feat / (2^b)
    as.integer(ceiling(a * exp(-k / b)))
  }


  # ------------------------------------------------------------------
  # .adaptive_reweighted_sampling
  # Samples `num_features` variable indices without replacement,
  # with probabilities proportional to absolute PLS coefficients.
  #
  # @param weights       Numeric vector of PLS regression coefficients.
  # @param num_features  Number of variables to select.
  # ------------------------------------------------------------------
  .adaptive_reweighted_sampling <- function(weights, num_features) {
    weights_abs  <- abs(weights)
    weights_norm <- weights_abs / sum(weights_abs)

    sample(
      seq_along(weights),
      size    = num_features,
      replace = FALSE,
      prob    = weights_norm
    )
  }


  # ------------------------------------------------------------------
  # .get_optimal_components
  # Determines a safe number of PLS components given data dimensions.
  #
  # @param n_samples       Number of available samples.
  # @param n_features      Number of features in the current subset.
  # @param max_components  Hard upper cap supplied by the user.
  # ------------------------------------------------------------------
  .get_optimal_components <- function(n_samples, n_features,
                                      max_components) {
    max_allowed <- min(n_samples - 1L, n_features, max_components)
    max(1L, as.integer(max_allowed))
  }


  # ------------------------------------------------------------------
  # .pls_rmsecv
  # Computes RMSECV via k-fold cross-validation using a PLS model.
  #
  # @param X_sel    Feature matrix (all samples × selected features).
  # @param y        Response vector.
  # @param ncomp    Number of PLS components to use.
  # @param n_folds  Number of CV folds.
  # @param seed     Optional integer seed for reproducible fold assignment.
  # ------------------------------------------------------------------
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
        preds     <- predict(model,
                             newdata = data.frame(X_test),
                             ncomp   = ncomp)[, 1L, 1L]
        sq_errors <- c(sq_errors, (y_test - preds)^2)
      }, error = function(e) NULL)
    }

    if (length(sq_errors) == 0L) return(Inf)
    sqrt(mean(sq_errors))
  }


  # ------------------------------------------------------------------
  # .run_monte_carlo
  # Executes all N Monte Carlo runs for one CARS iteration.
  # Returns the feature subset and RMSECV of the best run, or NULL
  # if every run fails.
  #
  # @param X                Full data matrix.
  # @param y                Response vector.
  # @param current_features Integer vector of active feature indices (1-based).
  # @param n_select         Target number of features for this iteration.
  # @param max_components   PLS component cap.
  # @param iteration        Current iteration index (used for seed offsets).
  # ------------------------------------------------------------------
  .run_monte_carlo <- function(X, y, current_features,
                               n_select, max_components, iteration) {

    n_samples          <- nrow(X)
    n_current_features <- length(current_features)
    iteration_rmsecv   <- numeric(0)
    iteration_features <- list()

    for (run in seq_len(N)) {

      # Sub-sample 80 % of data as calibration set
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
        coef_weights <- as.numeric(stats::coef(pls_model, ncomp = ncomp_cal))

        # Feature selection via Adaptive Reweighted Sampling
        selected_features <-
          if (n_current_features > n_select) {
            sel_idx <- .adaptive_reweighted_sampling(coef_weights, n_select)
            current_features[sel_idx]
          } else {
            current_features
          }

        # Evaluate subset via cross-validated RMSECV
        X_selected <- X[, selected_features, drop = FALSE]
        ncomp_cv   <- .get_optimal_components(n_samples,
                                              length(selected_features),
                                              max_components)

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
        message(sprintf("  [cars] Warning — Iter %d, Run %d: %s",
                        iteration, run, conditionMessage(e)))
      })
    }

    if (length(iteration_rmsecv) == 0L) return(NULL)

    best_idx <- which.min(iteration_rmsecv)
    list(
      features = iteration_features[[best_idx]],
      rmsecv   = iteration_rmsecv[[best_idx]]
    )
  }


  # ------------------------------------------------------------------ #
  #  PUBLIC method: fit()                                               #
  # ------------------------------------------------------------------ #

  #' @description
  #' Fit the CARS object to the high dimensional data `X` and response `y`.
  #'
  #' @param X              Numeric matrix of predictors (n_samples × n_features).
  #' @param y              Numeric response vector of length n_samples.
  #' @param max_components Integer cap on the number of PLS latent components Default `10`.
  #' @param save_plot      Boolean value, whether to generate plot Default `True`
  #' @param plot_path      Physical path to save the generated plot in local machine
  #'
  #' @return Named list with elements `best_features`, `best_rmsecv`,
  #'   `rmsecv_history`, and `num_features_history`.
  fit <- function(X,
                  y,
                  max_components = 10L,
                  plot = TRUE,
                  plot_path = "./cars_rmsecv_curve.jpg") {

    # ---- validate inputs ----
    if (!is.matrix(X))  X <- as.matrix(X)
    if (!is.numeric(X)) stop("`X` must be a numeric matrix.")
    if (!is.numeric(y)) stop("`y` must be a numeric vector.")
    if (nrow(X) != length(y))
      stop("`X` and `y` must have the same number of observations.")
    max_components <- as.integer(max_components)

    # Reset master seed to ensure reproducibility across repeated calls
    set.seed(random_state)

    n_features <- ncol(X)

    best_rmsecv          <- Inf
    best_features     <- NULL
    rmsecv_history       <- numeric(0)
    num_features_history <- integer(0)
    current_features     <- seq_len(n_features)

    message("[cars] Starting CARS — ", n_features, " initial features, ",
            max_iter, " max iterations, ", N, " MC runs each.")

    for (iteration in seq_len(max_iter)) {

      # Number of variables to keep this iteration (exponential schedule)
      n_select <- .exponential_decreasing_function(
        k      = iteration - 1L,   # 0-based index matches Python
        N_feat = n_features
      )

      if (n_select < 2L) {
        message(sprintf("[cars] Stopping at iteration %d: budget < 2 features.",
                        iteration))
        break
      }

      # Monte Carlo runs for this iteration
      mc_result <- .run_monte_carlo(
        X                = X,
        y                = y,
        current_features = current_features,
        n_select         = n_select,
        max_components   = max_components,
        iteration        = iteration
      )

      if (is.null(mc_result)) {
        message(sprintf("[cars] No successful MC runs at iteration %d. Stopping.",
                        iteration))
        break
      }

      rmsecv_history       <- c(rmsecv_history,       mc_result$rmsecv)
      num_features_history <- c(num_features_history, length(mc_result$features))

      if (mc_result$rmsecv < best_rmsecv) {
        best_rmsecv      <- mc_result$rmsecv
        best_features <- mc_result$features
      }

      current_features <- mc_result$features

      message(sprintf("[cars]  Iter %3d | features: %4d | RMSECV: %.4f",
                      iteration, length(current_features), mc_result$rmsecv))
    }

    if (!is.null(best_features)) {
      best_features <- sort(best_features)
    }

    # ── Plot ────────────────────────────────────────────────────────────────────
    best_idx        <- which.min(rmsecv_history)
    best_feature_pt <- num_features_history[best_idx]

    plot_df <- data.frame(
      num_features = num_features_history,
      rmsecv       = rmsecv_history
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = num_features, y = rmsecv)) +
      ggplot2::geom_line(color = "black", linewidth = 0.8) +
      ggplot2::geom_point(color = "black", size = 2) +
      ggplot2::geom_point(
        data = data.frame(x = best_feature_pt, y = best_rmsecv),
        ggplot2::aes(x = x, y = y),
        color = "darkred", shape = 15, size = 3,
        inherit.aes = FALSE
      ) +
      ggplot2::labs(
        x     = "Number of variables",
        y     = "RMSECV",
        title = paste0("Selected variables: ", length(best_wavelengths),
                       "  (RMSECV: ", round(best_rmsecv, 4), ")")
      ) +
      ggplot2::theme_minimal()

    print(p)

    if (save_plot) {
      dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
      ggplot2::ggsave(filename = plot_path, plot = p, width = 10, height = 6, dpi = 300)
      message("Plot saved to: ", plot_path)
    }

    structure(
      list(
        best_features        = best_features,
        best_rmsecv          = best_rmsecv,
        rmsecv_history       = rmsecv_history,
        num_features_history = num_features_history)
    )
  }

  # ---- return public interface only ----
  list(fit = fit)
}
