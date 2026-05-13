#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif

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
    
    // Fill Sigma upper triangle (row-major order)
    // numpy: Sigma[np.triu_indices(K)] yields row-major order:
    // Row 0: (0,0),(0,1),(0,2),...
    // Row 1: (1,1),(1,2),...
    int idx = 0;
    for (int r = 0; r < K; r++) {
      for (int c = r; c < K; c++) {
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

//' Compute parameter covariance upper triangles
//'
//' @param cov_vec Numeric vector containing concatenated upper triangles of covariance matrices.
//' @param offsets Integer vector of byte offsets (0-based) into the original file (needs conversion to index).
//' @param A The transformation matrix (n_params x K).
//' @param K Dimension of the covariance matrix (e.g. 17).
//' @return NumericMatrix of upper-triangle covariance (n_snps x n_cov_params).
// [[Rcpp::export]]
NumericMatrix compute_param_cov(NumericVector cov_vec, NumericVector offsets, NumericMatrix A, int K) {
  int n_snps = offsets.size();
  int n_params = A.nrow();
  int n_cov_elements = K * (K + 1) / 2;
  int n_cov_params = n_params * (n_params + 1) / 2;

  arma::mat A_mat(A.begin(), n_params, K, false);
  arma::mat Sigma(K, K, arma::fill::zeros);
  NumericMatrix out_cov(n_snps, n_cov_params);

  for (int i = 0; i < n_snps; i++) {
    double byte_offset = offsets[i];
    long long start_idx = (long long)(byte_offset / 4);

    if (start_idx < 0 || start_idx + n_cov_elements > cov_vec.size()) {
      for (int p = 0; p < n_cov_params; p++) out_cov(i, p) = NA_REAL;
      continue;
    }

    int idx = 0;
    for (int r = 0; r < K; r++) {
      for (int c = r; c < K; c++) {
        Sigma(r, c) = cov_vec[start_idx + idx];
        idx++;
      }
    }
    for (int r = 0; r < K; r++) {
      for (int c = r + 1; c < K; c++) {
        Sigma(c, r) = Sigma(r, c);
      }
    }

    arma::mat Vtheta = A_mat * Sigma * A_mat.t();

    int out_idx = 0;
    for (int r = 0; r < n_params; r++) {
      for (int c = r; c < n_params; c++) {
        out_cov(i, out_idx) = Vtheta(r, c);
        out_idx++;
      }
    }
  }

  return out_cov;
}

//' Compute GLS Stage 2 Parameter GWAS
//'
//' Uses the per-SNP tau covariance matrix from Stage 1 to estimate SNP-specific
//' generalized least-squares parameter effects and their covariance.
//'
//' @param beta_stage1 Tau-level SNP effects as a K x n_snps matrix.
//' @param cov_vec Numeric vector containing concatenated upper triangles of covariance matrices.
//' @param offsets Integer vector of byte offsets (0-based) into the original file.
//' @param W Weight/design matrix with K rows and n_params columns.
//' @param return_cov Whether to return parameter covariance upper triangles.
//' @param gls_ridge Relative diagonal ridge added only if the covariance solve fails.
//' @param n_threads Number of OpenMP threads.
//' @return List with params, se, Q, and optionally param_cov_upper.
// [[Rcpp::export]]
List compute_gls_param_gwas(const arma::mat& beta_stage1,
                            NumericVector cov_vec,
                            NumericVector offsets,
                            NumericMatrix W,
                            bool return_cov = false,
                            double gls_ridge = 1e-8,
                            int n_threads = 1) {
  const int K = beta_stage1.n_rows;
  const int n_snps = beta_stage1.n_cols;
  const int n_params = W.ncol();
  const int n_cov_elements = K * (K + 1) / 2;
  const int n_cov_params = n_params * (n_params + 1) / 2;

  if (W.nrow() != K) {
    stop("W must have the same number of rows as beta_stage1.");
  }
  if (offsets.size() < n_snps) {
    stop("offsets must have at least one entry per SNP.");
  }

  arma::mat W_mat(W.begin(), K, n_params, false);
  arma::mat params(n_params, n_snps, arma::fill::value(NA_REAL));
  arma::mat se(n_params, n_snps, arma::fill::value(NA_REAL));
  NumericVector Q(n_snps, NA_REAL);
  NumericMatrix param_cov_upper;
  if (return_cov) {
    param_cov_upper = NumericMatrix(n_snps, n_cov_params);
    std::fill(param_cov_upper.begin(), param_cov_upper.end(), NA_REAL);
  }

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (int i = 0; i < n_snps; i++) {
    const arma::vec beta_vec = beta_stage1.col(i);
    if (!beta_vec.is_finite()) continue;

    const long long start_idx = (long long)(offsets[i] / 4);
    if (start_idx < 0 || start_idx + n_cov_elements > cov_vec.size()) continue;

    arma::mat Sigma(K, K, arma::fill::zeros);
    bool cov_ok = true;
    int idx = 0;
    for (int r = 0; r < K; r++) {
      for (int c = r; c < K; c++) {
        double v = cov_vec[start_idx + idx];
        if (!std::isfinite(v)) {
          cov_ok = false;
          break;
        }
        Sigma(r, c) = v;
        Sigma(c, r) = v;
        idx++;
      }
      if (!cov_ok) break;
    }
    if (!cov_ok) continue;

    arma::mat rhs(K, n_params + 1);
    rhs.cols(0, n_params - 1) = W_mat;
    rhs.col(n_params) = beta_vec;

    arma::mat SinvRhs;
    bool ok = arma::solve(SinvRhs, Sigma, rhs,
                          arma::solve_opts::fast + arma::solve_opts::likely_sympd);

    if ((!ok || !SinvRhs.is_finite()) && gls_ridge > 0.0) {
      const arma::vec d = Sigma.diag();
      double diag_scale = arma::mean(d.elem(arma::find_finite(d)));
      if (!std::isfinite(diag_scale) || diag_scale <= 0.0) diag_scale = 1.0;
      arma::mat Sigma_ridge = Sigma;
      Sigma_ridge.diag() += gls_ridge * diag_scale;
      ok = arma::solve(SinvRhs, Sigma_ridge, rhs,
                       arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      if (ok && SinvRhs.is_finite()) {
        Sigma = Sigma_ridge;
      }
    }

    if (!ok || !SinvRhs.is_finite()) continue;

    arma::mat SinvW = SinvRhs.cols(0, n_params - 1);
    arma::vec SinvBeta = SinvRhs.col(n_params);
    arma::mat WtSinvW = W_mat.t() * SinvW;
    arma::vec WtSinvBeta = W_mat.t() * SinvBeta;

    arma::vec theta;
    ok = arma::solve(theta, WtSinvW, WtSinvBeta,
                     arma::solve_opts::fast + arma::solve_opts::likely_sympd);
    if (!ok || !theta.is_finite()) continue;

    arma::mat Vtheta;
    ok = arma::solve(Vtheta, WtSinvW, arma::eye(n_params, n_params),
                     arma::solve_opts::fast + arma::solve_opts::likely_sympd);
    if (!ok || !Vtheta.is_finite()) continue;

    arma::vec resid = beta_vec - W_mat * theta;
    arma::vec SinvResid = SinvBeta - SinvW * theta;
    double q_val = arma::dot(resid, SinvResid);

    params.col(i) = theta;
    for (int p = 0; p < n_params; p++) {
      double var = Vtheta(p, p);
      se(p, i) = (std::isfinite(var) && var > 0.0) ? std::sqrt(var) : NA_REAL;
    }
    Q[i] = std::isfinite(q_val) ? q_val : NA_REAL;

    if (return_cov) {
      int out_idx = 0;
      for (int r = 0; r < n_params; r++) {
        for (int c = r; c < n_params; c++) {
          param_cov_upper(i, out_idx) = Vtheta(r, c);
          out_idx++;
        }
      }
    }
  }

  List out = List::create(
    _["params"] = params,
    _["se"] = se,
    _["Q"] = Q
  );
  if (return_cov) {
    out["param_cov_upper"] = param_cov_upper;
  }
  return out;
}
