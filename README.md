
# Competitive Adaptive Reweighted Sampling (CARS)


<!-- badges: start -->
<!-- badges: end -->

> Competitive Adaptive Reweighted Sampling (CARS) Algorithm for variable selection in high dimensional dataset

## Installation

You can install the development version of cars like so:

``` r
install.packages(c("pls", "readxl", "ggplot2"))

# Install cars from local source
install.packages("path/to/cars_0.1.0.tar.gz", repos = NULL, type = "source")
or 
remotes::install_github("mah-iasri/cars")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(cars)
set.seed(1)
X <- matrix(rnorm(100 * 200), nrow = 100)
y <- X[, 5] * 2 + X[, 50] * -1.5 + rnorm(100, sd = 0.5)

# Instantiate and fit
cars <- CARSAlgorithm(max_iter = 15, N = 30, cv_folds = 5, random_state = 42)
result <- cars$fit(X, y, max_components = 8)

cat("Best RMSECV      :", result$best_rmsecv, "\n")
cat("Selected features:", result$best_features, "\n")
```

