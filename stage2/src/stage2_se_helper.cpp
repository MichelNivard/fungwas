#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Compute Calibrated Standard Errors
//' 
//' @param cov_vec Numeric vector containing concatenated upper triangles of covariance matrices.
//' @param offsets Integer vector of byte offsets (0-based) into the original file (needs conversion to index).
//' @param A The transformation matrix (n_params x K).
//' @param K Dimension of the covariance matrix (e.g. 17).
//' @return NumericMatrix of standard errors (n_snps x n_params).
// [[Rcpp::export]]
NumericMatrix compute_calibrated_se(NumericVector cov_vec, NumericVector offsets, NumericMatrix A, int K) {
  int n_snps = offsets.size();
  int n_params = A.nrow();
  int n_cov_elements = K * (K + 1) / 2;
  
  arma::mat A_mat(A.begin(), n_params, K, false); // Wrap existing memory
  arma::mat Sigma(K, K, arma::fill::zeros);
  NumericMatrix out_se(n_snps, n_params);
  
  // R uses 32-bit integers for vector indexing, but offsets might be large if file is >2GB.
  // However, cov_vec is in memory, so it fits in R vector (max 2^31-1 elements).
  // The offsets passed here are from the .idx file, which are BYTE offsets.
  // cov_vec is a vector of floats (4 bytes).
  
  for (int i = 0; i < n_snps; i++) {
    double byte_offset = offsets[i];
    
    // Convert byte offset to vector index (0-based)
    // dividing by 4 because floats are 4 bytes
    // R is 1-based, but C++ is 0-based.
    // offsets[i] is start of record in bytes.
    long long start_idx = (long long)(byte_offset / 4);
    
    if (start_idx < 0 || start_idx + n_cov_elements > cov_vec.size()) {
       // Index out of bounds, skip or produce NA
       for (int p=0; p<n_params; p++) out_se(i, p) = NA_REAL;
       continue;
    }
    
    // Fill Sigma upper triangle
    // R stores matrices in column-major order.
    // upper.tri(diag=TRUE) extraction follows this order:
    // Col 0: (0,0)
    // Col 1: (0,1), (1,1)
    // Col 2: (0,2), (1,2), (2,2)
    // ...
    
    int idx = 0;
    for (int c = 0; c < K; c++) {
      for (int r = 0; r <= c; r++) {
        Sigma(r, c) = cov_vec[start_idx + idx];
        idx++;
      }
    }
    // Symmetrize
    for (int r = 0; r < K; r++) {
      for (int c = r + 1; c < K; c++) {
        Sigma(c, r) = Sigma(r, c);
      }
    }
    
    // V_param = A * Sigma * A.t()
    // We only need the diagonal (variances)
    // Optimization: V_j = (A.row(j)) * Sigma * (A.row(j)).t()
    
    for (int p = 0; p < n_params; p++) {
       arma::rowvec a_row = A_mat.row(p);
       double var = arma::as_scalar(a_row * Sigma * a_row.t());
       if (var > 0) {
         out_se(i, p) = std::sqrt(var);
       } else {
         out_se(i, p) = 0.0;
       }
    }
  }
  
  return out_se;
}
