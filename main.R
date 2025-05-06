
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
        
        results <- run_simulation_grid(nrep = 50)
        print(results)
}

# Run the main function
main()
