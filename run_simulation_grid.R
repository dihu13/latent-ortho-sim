# run_simulation_grid.R
run_simulation_grid <- function(nrep = 50,
                                results_path = "results",
                                d0 = c(25, 100),
                                p0 = c(1, 5),
                                n0 = c(200, 500),
                                q0 = c(2, 5),
                                pctB0 = c(0.1, 0.9)) {
        dir.create(results_path, showWarnings = FALSE)
        grid <- expand.grid(n = n0, d = d0, p = p0, q = q0, pct.B = pctB0)
        mse_summary <- data.frame()
        
        for (i in seq_len(nrow(grid))) {
                n <- grid$n[i]; d <- grid$d[i]; p <- grid$p[i]; q <- grid$q[i]; pct.B <- grid$pct.B[i]
                
                mse_naive <- mse_oracle <- mse_adj <- coverage_adj <- numeric(nrep)
                
                for (j in seq_len(nrep)) {
                        dat <- simulate_data(n, p, q, d, pct.B)
                        B_true <- dat$B
                        
                        B_naive <- estimate_naive(dat$X, dat$Y)
                        B_oracle <- estimate_oracle(dat$X, dat$Z, dat$Y)
                        est_adj <- estimate_factor_adjusted(dat$X, dat$Y, q, d)
                        
                        if (is.null(est_adj)) next
                        
                        mse_naive[j] <- mean((B_naive - B_true)^2)
                        mse_oracle[j] <- mean((B_oracle - B_true)^2)
                        mse_adj[j] <- mean((est_adj$B_hat - B_true)^2)
                        
                        if (!is.null(est_adj$seB)) {
                                z_val <- qnorm(0.975)
                                lower <- est_adj$B_hat - z_val * est_adj$seB
                                upper <- est_adj$B_hat + z_val * est_adj$seB
                                coverage <- mean(B_true >= lower & B_true <= upper)
                                coverage_adj[j] <- coverage
                        }
                }
                
                res <- data.frame(n = n, d = d, p = p, q = q, pctB = pct.B,
                                  mse_naive = mean(mse_naive, na.rm = TRUE),
                                  mse_oracle = mean(mse_oracle, na.rm = TRUE),
                                  mse_adj = mean(mse_adj, na.rm = TRUE),
                                  coverage_adj = mean(coverage_adj, na.rm = TRUE))
                
                fname <- paste0("res_", n, "_", d, "_", p, "_", q, "_", pct.B, ".rds")
                saveRDS(res, file.path(results_path, fname))
                
                mse_summary <- rbind(mse_summary, res)
        }
        write.csv(mse_summary, file.path(results_path, "summary.csv"), row.names = FALSE)
        return(mse_summary)
}

