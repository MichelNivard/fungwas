#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix fit_multi_models(const arma::mat& beta_stage1,
                               const arma::mat& se_stage1,
                               const List& W_list,
                               double N,
                               int n_threads = 1) {
  const uword P = beta_stage1.n_rows;
  const uword T = beta_stage1.n_cols;
  const int n_models = W_list.size();

  if (se_stage1.n_rows != P || se_stage1.n_cols != T) {
    stop("se_stage1 must have the same shape as beta_stage1.");
  }
  if (n_models == 0) {
    return NumericMatrix(P, 0);
  }

  std::vector<arma::mat> W_mats;
  std::vector<uword> k_sizes;
  W_mats.reserve(n_models);
  k_sizes.reserve(n_models);

  for (int m = 0; m < n_models; m++) {
    arma::mat W = as<arma::mat>(W_list[m]);
    if (W.n_rows != T) {
      stop("Each W matrix must have T rows.");
    }
    W_mats.push_back(W);
    k_sizes.push_back(W.n_cols);
  }

  NumericMatrix out(P, 2 * n_models);

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (uword i = 0; i < P; i++) {
    arma::rowvec beta = beta_stage1.row(i);
    arma::rowvec se = se_stage1.row(i);
    arma::rowvec weight = 1.0 / arma::square(se);

    for (uword t = 0; t < T; t++) {
      if (!std::isfinite(weight(t)) || !std::isfinite(beta(t))) {
        weight(t) = 0.0;
      }
    }

    arma::vec beta_vec = beta.t();
    arma::vec weight_vec = weight.t();

    for (int m = 0; m < n_models; m++) {
      const arma::mat& W = W_mats[m];
      uword k = k_sizes[m];

      arma::mat W_weighted = W.each_col() % weight_vec;
      arma::mat WtSinvW = W.t() * W_weighted;
      arma::vec WtSinvBeta = W.t() * (weight_vec % beta_vec);

      arma::vec theta;
      bool ok = arma::solve(theta, WtSinvW, WtSinvBeta,
                            arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      if (!ok || !theta.is_finite()) {
        out(i, 2 * m) = NA_REAL;
        out(i, 2 * m + 1) = NA_REAL;
        continue;
      }

      arma::vec resid = beta_vec - W * theta;
      double Q = arma::dot(resid % weight_vec, resid);
      out(i, 2 * m) = Q + 2.0 * k;
      out(i, 2 * m + 1) = Q + k * std::log(N);
    }
  }

  return out;
}
