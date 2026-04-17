
# Competitive Adaptive Reweighted Sampling (CARS)


<!-- badges: start -->
<!-- badges: end -->

> Competitive Adaptive Reweighted Sampling (CARS) Algorithm

## Description

The carsAlgo implements Competitive Adaptive Reweighted Sampling (CARS) 
    algorithm for variable selection from high-dimensional dataset using in Partial 
    Least Squares (PLS) regression models. CARS algorithm iteratively applies the 
    Monte Carlo sub-sampling and exponential variable elimination techniques to 
    identify/select the most informative variables/features subjected to minimal 
    cross-validated RMSE score. The implementation of CARS algorithm is inspired 
    from the work of Li et al. (2009) <doi:10.1016/j.aca.2009.06.046>. 
    This algorithm is widely applied in near-infrared (NIR), mid-infrared (MIR), 
    hyperspectral chemometrics areas, etc.

## Installation

You can install the development version of cars like so:

``` r
install.packages(c("pls", "ggplot2", "stats"))

# Install cars from local source
install.packages("path/to/carsAlgo_0.5.0.tar.gz", repos = NULL, type = "source")
or 
remotes::install_github("mah-iasri/carsAlgo")
```

## Example

This is a basic example of how cars can be used:

``` r
library(carsAlgo)
set.seed(1)
X <- matrix(rnorm(100 * 200), nrow = 100)
y <- X[, 5] * 2 + X[, 50] * -1.5 + rnorm(100, sd = 0.5)

cars_obj <- CARSAlgorithm(max_iter = 15, N = 30, cv_folds = 5)
result   <- fit(cars_obj, X, y, max_components = 8,
                plot_path = file.path(tempdir(), "cars_rmsecv_curve.jpg"))

cat("Best RMSECV      :", result$best_rmsecv, "\n")
cat("Selected features:", result$best_features, "\n")
```

## Reference

Li, H., Liang, Y., Xu, Q., & Cao, D. (2009). Key wavelengths screening using
competitive adaptive reweighted sampling method for multivariate calibration.
*Analytica Chimica Acta*, 648(1), 77–84. <https://doi.org/10.1016/j.aca.2009.06.046>

## Citation

Haque, M.A., Ghosh, A., Karamakar, S., Sachan, H. and Kumari, S. (2026). carsAlgo. R Package. *CRAN*. Doi: 10.32614/CRAN.package.carsAlgo


