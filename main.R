
# main.R
main <- function() {
        library(MASS)
        library(Matrix)
        library(statmod)
        library(Rcpp)
        library(RcppArmadillo)
        
        source("simulate_data.R")
        source("estimate_methods.R")
        source("run_simulation_grid.R")
        
        # Placeholder for Rcpp function
        rcpp_cov_fn <- NULL
        if (file.exists("ortho_cov_central.so")) {
                dyn.load("ortho_cov_central.so")
                rcpp_cov_fn <- get("CovarianceCpp_ortho1")
        }
        
        results <- run_simulation_grid(nrep = 20)
        print(head(results))
}

# Run the main function
main()
