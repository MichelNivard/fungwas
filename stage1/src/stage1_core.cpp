// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List openmp_info() {
  #ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  #else
  int max_threads = 1;
  #endif
  return List::create(
    Named("openmp_enabled") = LogicalVector::create(
      #ifdef _OPENMP
      true
      #else
      false
      #endif
    ),
    Named("max_threads") = max_threads
  );
}

// [[Rcpp::export]]
List stage1_rif_ols(const arma::mat& G,
                    const arma::mat& RIF_resid,
                    const arma::mat& Q,
                    double min_mac = 5.0,
                    double min_var = 1e-8,
                    int n_threads = 1) {
  const uword N = G.n_rows;
  const uword M = G.n_cols;
  const uword T = RIF_resid.n_cols;

  mat betas(M, T, fill::value(datum::nan));
  mat se_hc0(M, T, fill::value(datum::nan));

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (uword j = 0; j < M; j++) {
    vec g = G.col(j);

    uword n_valid = 0;
    double g_sum = 0.0;
    for (uword i = 0; i < N; i++) {
      if (is_finite(g(i))) {
        g_sum += g(i);
        n_valid++;
      }
    }

    if (n_valid < 100) continue;

    double g_mean = g_sum / n_valid;
    double mac = std::min(g_sum, 2.0 * n_valid - g_sum);
    if (mac < min_mac) continue;

    for (uword i = 0; i < N; i++) {
      if (!is_finite(g(i))) {
        g(i) = g_mean;
      }
    }

    vec Qtg = Q.t() * g;
    vec g_resid = g - Q * Qtg;

    double g2 = dot(g_resid, g_resid);
    if (g2 < min_var) continue;

    vec g_sq = g_resid % g_resid;

    for (uword t = 0; t < T; t++) {
      double g_rif = dot(g_resid, RIF_resid.col(t));
      double beta = g_rif / g2;
      betas(j, t) = beta;

      double hc0_sum = 0.0;
      const double* rif_ptr = RIF_resid.colptr(t);
      for (uword i = 0; i < N; i++) {
        double e_i = rif_ptr[i] - g_resid(i) * beta;
        hc0_sum += g_sq(i) * e_i * e_i;
      }
      se_hc0(j, t) = std::sqrt(hc0_sum) / g2;
    }
  }

  return List::create(
    Named("betas") = betas,
    Named("se_hc0") = se_hc0
  );
}

// [[Rcpp::export]]
List stage1_rif_ols_multi(const arma::mat& G,
                          const List& rif_list,
                          const arma::mat& Q,
                          double min_mac = 5.0,
                          double min_var = 1e-8,
                          int n_threads = 1) {
  const uword N = G.n_rows;
  const uword M = G.n_cols;
  const int n_pheno = rif_list.size();

  if (n_pheno == 0) {
    return List::create();
  }

  std::vector<NumericMatrix> rif_inputs;
  rif_inputs.reserve(n_pheno);
  std::vector<arma::mat> rif_mats;
  rif_mats.reserve(n_pheno);
  std::vector<arma::mat> betas_list;
  std::vector<arma::mat> se_list;
  betas_list.reserve(n_pheno);
  se_list.reserve(n_pheno);

  for (int p = 0; p < n_pheno; p++) {
    NumericMatrix rif = rif_list[p];
    if (rif.nrow() != static_cast<int>(N)) {
      stop("All RIF matrices must have the same number of rows as G.");
    }
    rif_inputs.push_back(rif);
    rif_mats.emplace_back(rif.begin(), rif.nrow(), rif.ncol(), false, true);
    betas_list.emplace_back(M, rif.ncol(), fill::value(datum::nan));
    se_list.emplace_back(M, rif.ncol(), fill::value(datum::nan));
  }

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (uword j = 0; j < M; j++) {
    vec g = G.col(j);

    uword n_valid = 0;
    double g_sum = 0.0;
    for (uword i = 0; i < N; i++) {
      if (is_finite(g(i))) {
        g_sum += g(i);
        n_valid++;
      }
    }

    if (n_valid < 100) continue;

    double g_mean = g_sum / n_valid;
    double mac = std::min(g_sum, 2.0 * n_valid - g_sum);
    if (mac < min_mac) continue;

    for (uword i = 0; i < N; i++) {
      if (!is_finite(g(i))) {
        g(i) = g_mean;
      }
    }

    vec Qtg = Q.t() * g;
    vec g_resid = g - Q * Qtg;

    double g2 = dot(g_resid, g_resid);
    if (g2 < min_var) continue;

    vec g_sq = g_resid % g_resid;

    for (int p = 0; p < n_pheno; p++) {
      const arma::mat& rif = rif_mats[p];
      arma::mat& betas = betas_list[p];
      arma::mat& se_hc0 = se_list[p];
      const uword T = rif.n_cols;

      for (uword t = 0; t < T; t++) {
        double g_rif = dot(g_resid, rif.col(t));
        double beta = g_rif / g2;
        betas(j, t) = beta;

        double hc0_sum = 0.0;
        const double* rif_ptr = rif.colptr(t);
        for (uword i = 0; i < N; i++) {
          double e_i = rif_ptr[i] - g_resid(i) * beta;
          hc0_sum += g_sq(i) * e_i * e_i;
        }
        se_hc0(j, t) = std::sqrt(hc0_sum) / g2;
      }
    }
  }

  List results(n_pheno);
  for (int p = 0; p < n_pheno; p++) {
    results[p] = List::create(
      Named("betas") = betas_list[p],
      Named("se_hc0") = se_list[p]
    );
  }

  return results;
}

// [[Rcpp::export]]
List stage1_block_scores(const arma::mat& G,
                         const arma::mat& RIF_resid,
                         const arma::mat& Q,
                         const arma::uvec& block_ids,
                         int n_blocks,
                         double min_mac = 5.0,
                         double min_var = 1e-8,
                         int n_threads = 1) {
  const uword N = G.n_rows;
  const uword M = G.n_cols;
  const uword T = RIF_resid.n_cols;

  // Output: M SNPs, each with a matrix of (n_blocks x (1 + T))
  // Col 0: sum(G^2), Cols 1-T: sum(G*Y)
  field<mat> block_stats(M);

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (uword j = 0; j < M; j++) {
    vec g = G.col(j);
    uword n_valid = 0;
    double g_sum = 0.0;
    for (uword i = 0; i < N; i++) {
      if (is_finite(g(i))) {
        g_sum += g(i);
        n_valid++;
      }
    }

    if (n_valid < 100) continue;

    double g_mean = g_sum / n_valid;
    for (uword i = 0; i < N; i++) {
      if (!is_finite(g(i))) g(i) = g_mean;
    }

    vec Qtg = Q.t() * g;
    vec g_resid = g - Q * Qtg;

    mat stats(n_blocks, 1 + T, fill::zeros);

    for (uword i = 0; i < N; i++) {
      uword b = block_ids(i);
      if (b >= n_blocks) continue; // safety

      double g_i = g_resid(i);
      stats(b, 0) += g_i * g_i;
      for (uword t = 0; t < T; t++) {
        stats(b, 1 + t) += g_i * RIF_resid(i, t);
      }
    }
    block_stats(j) = stats;
  }

  return List::create(Named("block_stats") = block_stats);
}

// [[Rcpp::export]]
List ols_block_scores(const arma::mat& G,
                      const arma::vec& Y_resid,
                      const arma::mat& Q,
                      const arma::uvec& block_ids,
                      int n_blocks,
                      double min_mac = 5.0,
                      double min_var = 1e-8,
                      int n_threads = 1) {
  // OLS version: single outcome instead of multiple taus
  // Returns block-level sums for jackknife variance estimation
  
  const uword N = G.n_rows;
  const uword M = G.n_cols;

  // Output matrices: M SNPs x n_blocks
  // D_block[j,b] = sum of G_resid^2 in block b for SNP j
  // N_block[j,b] = sum of G_resid * Y_resid in block b for SNP j
  mat D_block(M, n_blocks, fill::zeros);
  mat N_block(M, n_blocks, fill::zeros);

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (uword j = 0; j < M; j++) {
    vec g = G.col(j);
    
    // Count valid and compute mean
    uword n_valid = 0;
    double g_sum = 0.0;
    for (uword i = 0; i < N; i++) {
      if (is_finite(g(i))) {
        g_sum += g(i);
        n_valid++;
      }
    }

    if (n_valid < 100) continue;

    // Mean impute missing
    double g_mean = g_sum / n_valid;
    for (uword i = 0; i < N; i++) {
      if (!is_finite(g(i))) g(i) = g_mean;
    }

    // Residualize genotype on covariates
    vec Qtg = Q.t() * g;
    vec g_resid = g - Q * Qtg;

    // Accumulate block-level statistics
    for (uword i = 0; i < N; i++) {
      uword b = block_ids(i);
      if (b >= (uword)n_blocks) continue;

      double g_i = g_resid(i);
      D_block(j, b) += g_i * g_i;
      N_block(j, b) += g_i * Y_resid(i);
    }
  }

  return List::create(
    Named("D_block") = D_block,
    Named("N_block") = N_block
  );
}
