
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Enable OpenMP
#define ARMA_USE_OPENMP
#include <RcppArmadillo.h>

// Enable LAPACK
#define ARMA_USE_LAPACK
#include <RcppArmadillo.h>

// Enable BLAS
#define ARMA_USE_BLAS
#include <RcppArmadillo.h>



// [[Rcpp::export]]
List rcpp_cov_fn(arma::mat X,
                          arma::mat Y,
                          arma::mat D_hat,
                          arma::vec Psi_hat,
                          int ncores){
        int n = X.n_rows;
        int p = X.n_cols;
        int q = D_hat.n_rows;
        int d = D_hat.n_cols;
        
        arma::mat MX = arma::inv_sympd(X.t() * X);
        arma::mat Sigma_tilt_Y_hat = arma::diagmat(Psi_hat) + arma::trans(D_hat) * D_hat;
        arma::mat hat_CovBs = arma::kron(Sigma_tilt_Y_hat, MX);
        arma::mat ImMJ = arma::eye<arma::mat>(d, d) -  arma::ones(d, d) / d;
        
        arma::mat inv_Sigma_tilt_Y_hat =arma::inv_sympd(Sigma_tilt_Y_hat);
        arma::mat DSD = D_hat * inv_Sigma_tilt_Y_hat * D_hat.t();
        arma::mat SD = inv_Sigma_tilt_Y_hat * D_hat.t();
        
        // Covariance matrix of D
        
        arma::mat I_dd(d * q, d * q);
        arma::mat I_dp(d * q, d);
        arma::mat I_pp(d, d);
        
        int size = d*q;
        
#if defined(_OPENMP)
#pragma omp parallel num_threads(ncores)
#pragma omp for
#endif
        
        for (int u = 0; u < size; u++) {
                for (int v = u; v < size; v++) {
                        int i = (u / q) + 1;
                        int r = u - (i - 1) * q + 1;
                        int j = (v / q) + 1;
                        int s = v - (j - 1) * q + 1;
                        
                        I_dd(u, v) = inv_Sigma_tilt_Y_hat(i - 1, j - 1) * DSD(r - 1, s - 1) +
                                SD(i - 1, s - 1) * SD(j - 1, r - 1);
                        // Filling the lower triangle as well
                        I_dd(v, u) = I_dd(u, v);
                }
        }
        for (int i = 0; i < d; i++) {
                for (int r = 0; r < q; r++) {
                        for (int j = 0; j < d; j++) {
                                I_dp(i * q + r, j) = inv_Sigma_tilt_Y_hat(i, j) * SD(j, r);
                        }
                }
        }
        
        for (int i = 0; i < d; i++) {
                for (int j = i; j < d; j++) {
                        I_pp(i, j) = 0.5 * inv_Sigma_tilt_Y_hat(i, j) * inv_Sigma_tilt_Y_hat(i, j);
                        // Filling the lower triangle as well
                        I_pp(j, i) = I_pp(i, j);
                }
        }
        
        
        arma::mat a_gd(q * d, q * (q - 1) / 2);
        // Compute a_gd
        for (int u = 0; u < (q - 1); u++) {
                for (int v = u + 1; v < q; v++) {
                        for (int i = 0; i < d; i++) {
                                for (int r = 0; r < q; r++) {
                                        
                                        int row_index = i * q + r;
                                        int col_index = (u * (2 * q - u - 1) / 2) + (v - u - 1);
                                        a_gd(row_index, col_index) = ((r == u) * D_hat(v, i) + (r == v) * D_hat(u, i)) / Psi_hat(i);
                                }
                        }
                }
        }
        
        
        arma::mat a_gp(d, q * (q - 1) / 2);
        // Compute a_gp
        for (int u = 0; u < (q - 1); u++) {
                for (int v = u + 1; v < q; v++) {
                        for (int i = 0; i < d; i++) {
                                int row_index = i;
                                int col_index = (u * (2 * q - u - 1) / 2) + (v - u - 1);
                                a_gp(row_index, col_index) = -D_hat(u, i) * D_hat(v, i) / (Psi_hat(i) * Psi_hat(i));
                        }
                }
        }
        
        
        
        // Calculate the information matrix
        
        // join columns (=vertically, i.e. matrices must have the same number of cols)
        arma::mat r1 = arma::join_rows(I_dd, I_dp, a_gd);
        arma::mat r2 = arma::join_rows(I_dp.t(),I_pp,a_gp);
        arma::mat z0 = arma::zeros(q * (q - 1) / 2, q * (q - 1) / 2);
        arma::mat r3 = arma::join_rows(a_gd.t(),a_gp.t(),z0);  
        
        //Rcpp::Rcout << "r1 dimensions: " << r1.n_rows << " x " << r1.n_cols << std::endl;
        //Rcpp::Rcout << "r2 dimensions: " << r2.n_rows << " x " << r2.n_cols << std::endl;
        //Rcpp::Rcout << "r3 dimensions: " << r3.n_rows << " x " << r3.n_cols << std::endl;
        
        arma::mat I  = arma::join_cols( r1, r2, r3 );
        
        // The estimate of covariance matrix for D, Psi
        
        
        arma::mat n_hat_Cov = inv(I);
        arma::mat hat_Cov = n_hat_Cov/n;
        arma::mat hat_CovD = hat_Cov.submat(0,0,q*d-1,q*d-1);
        arma::mat hat_CovPsi = hat_Cov.submat(q*d,q*d,q*d+d-1,q*d+d-1);
        arma::mat hat_seD = arma::reshape(arma::sqrt(hat_CovD.diag()),q,d);
        
        // Orthogonal Scenario
        arma::mat Q = D_hat*ImMJ*D_hat.t();
        arma::mat QD = arma::inv_sympd(Q)*D_hat;
        arma::mat DQD = D_hat.t()*QD;
        arma::mat ImHD = arma::eye<arma::mat>(d, d) - DQD*ImMJ;
        arma::mat B2 = MX*X.t()*Y;
        
        // The gradient matrix of D over B
        
        arma::mat aBD(p * d, q * d);
        for (int i = 0; i < q; i++) {
                for (int j = 0; j < d; j++) {
                        arma::mat aDd = arma::zeros<arma::mat>(q, d);
                        aDd(i, j) = 1;
                        arma::mat aBd = -B2 * ImMJ * (ImHD * aDd.t() * QD + arma::trans(ImHD * aDd.t() * QD));
                        aBD.col((j * q) + i) = arma::vectorise(aBd);
                }
        }
        // The estimate for covariance of B from D
        arma::mat hat_CovB_D = aBD * hat_CovD * aBD.t();
        
        // The gradient matrix of B over B_star
        arma::mat aBBs(p * d, p * d);
        for (int i = 0; i < p; i++) {
                for (int j = 0; j < d; j++) {
                        arma::mat aBb = arma::zeros<arma::mat>(p, d);
                        aBb(i, j) = 1;
                        arma::mat aBbs = aBb * ImHD.t();
                        aBBs.col((j * p) + i) = arma::vectorise(aBbs);
                }
        }
        
        // The estimate for covariance of B from D
        arma::mat hat_CovB_Bs = aBBs * hat_CovBs * aBBs.t();
        // The estimate of cov(B)
        arma::mat hat_CovB =  hat_CovB_D + hat_CovB_Bs;
        
        // Calculate the rate of true value dropping out of the CI
        arma::vec hat_VarB = arma::diagvec(hat_CovB);
        arma::mat hat_seB = arma::reshape(arma::sqrt(hat_VarB),p,d);
        
        // Create a List with results
        List result = List::create(
                Named("CovD") = hat_CovD,
                Named("seD") = hat_seD,
                Named("CovB") = hat_CovB,
                Named("seB") = hat_seB
        );
        return result;
        
}