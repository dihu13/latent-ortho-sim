# simulate_data.R
simulate_data <- function(n, p, q, d, pct.B) {
        require(MASS)
        require(invgamma)
        
        # Simulate covariates X
        if (p == 1) {
                X <- matrix(rnorm(n), ncol = 1)
        } else {
                cor0 <- matrix(runif(p, -0.7, 0.7), ncol = p)
                corrX <- t(cor0) %*% cor0
                diag(corrX) <- 0
                Sigma <- diag(1, p) + corrX
                X <- mvrnorm(n, rep(0, p), Sigma)
                X[, 1:3] <- (X[, 1:3] > 0) + 0
                X <- scale(X)
        }
        
        A <- matrix(round(rnorm(p * q), 1), ncol = q)
        W <- round(mvrnorm(n, rep(0, q), diag(1, q)), 1)
        B <- matrix(sample(c(2, -2), p * d, replace = TRUE), ncol = d)
        Psi <- rinvgamma(d, shape = 10, rate = 9)
        sampB <- sample(c(rep(1, ceiling(pct.B * d * p)), rep(0, floor((1 - pct.B) * d * p))))
        B[sampB == 0] <- 0
        E <- sapply(Psi, function(psi) rnorm(n, 0, sqrt(psi)))
        
        MJ <- matrix(1, nrow = d, ncol = d) / d
        ImMJ <- diag(1, d) - MJ
        base.orth.B <- svd(B %*% ImMJ, nv = ncol(B))$v[, (p + 1):d]
        M <- matrix(round(rnorm(q * (d - p)), 1), ncol = d - p)
        D0 <- t(M %*% t(base.orth.B))
        U <- eigen((t(D0 / Psi) %*% D0) / d)$vectors
        D <- t(U) %*% t(D0)
        
        Z <- X %*% A + W
        Y <- X %*% B + Z %*% D + E
        B_star <- B + A %*% D
        
        list(X = X, Y = Y, B = B, D = D, A = A, Z = Z, Psi = Psi, B_star = B_star)
}


