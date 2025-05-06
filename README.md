# Latent Confounder Simulation Study

This project simulates multivariate regression scenarios where the outcome is influenced by latent (unobserved) confounders. It evaluates several estimators for the primary effect of interest under varying conditions, including:

- **Naive OLS estimator** (ignores confounding),
- **Oracle estimator** (assume true latent confounders are observable),
- **Factor-adjusted estimator** (recovers latent structure via factor analysis and corrects for bias).

##  Features

- Modular, clean R codebase
- Grid-based simulation framework across combinations of sample size (`n`), predictors (`p`), outcomes (`d`), and latent dimension (`q`)
- Optional integration with high-performance C++ (via `RcppArmadillo`) for computing standard errors and confidence interval coverage
- Results saved per-scenario and summarized in a final `.csv`

##  File Structure

- `simulate_data.R`: Generates synthetic datasets with latent confounding
- `estimate_methods.R`: Implements three estimators
- `run_simulation_grid.R`: Runs simulations over user-specified parameter grid
- `main.R`: Entry point to run everything
- `ortho_cov.cpp`: Optional Rcpp file for efficient standard error estimation

## Requirements


```r
install.packages(c("MASS", "psych", "statmod", "Rcpp", "RcppArmadillo","cate", "invgamma"))