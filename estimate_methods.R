# estimate_methods.R
estimate_naive <- function(X, Y) {
        solve(t(X) %*% X) %*% t(X) %*% Y
}

estimate_oracle <- function(X, Z, Y) {
        P <- cbind(X, Z)
        model <- lm(Y ~ P - 1)
        coef(model)[1:ncol(X), ]  # return only coefficients for X
}

estimate_factor_adjusted <- function(X, Y, q, d) {
        require(cate)
        library(Rcpp)
        library(RcppArmadillo)
        sourceCpp("ortho_cov.cpp")
        B_init <- estimate_naive(X, Y)
        tilde_Y <- Y - X %*% B_init
        res <- try(factor.analysis(tilde_Y, r = q, method = "ml"), silent = TRUE)
        if (inherits(res, "try-error")) return(NULL)
        
        Gamma <- res$Gamma
        Psi_hat <- res$Sigma
        D_cate0 <- t(Gamma)
        MJ <- matrix(1, nrow = d, ncol = d) / d
        ImMJ <- diag(1, d) - MJ
        D_cate <- t(eigen((D_cate0 %*% (t(D_cate0) / Psi_hat)) / d)$vectors) %*% D_cate0
        
        A_hat <- B_init %*% ImMJ %*% t(D_cate) %*% solve(D_cate %*% ImMJ %*% t(D_cate))
        B_hat <- B_init - A_hat %*% D_cate
        
        if (!is.null(rcpp_cov_fn)) {
                se_result <- rcpp_cov_fn(X, Y, D_cate, Psi_hat, parallel::detectCores())
                seB <- se_result$seB
        } else {
                seB <- matrix(NA, nrow = ncol(X), ncol = d)  # placeholder
        }
        
        list(B_hat = B_hat, A_hat = A_hat, D_hat = D_cate, Psi_hat = Psi_hat, seB = seB)
} 