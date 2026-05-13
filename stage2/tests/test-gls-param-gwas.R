# ============================================================
# Test: tau_cov param_gwas uses GLS betas/SEs when covariance is available
# ============================================================

suppressPackageStartupMessages(library(Rcpp))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg)
} else {
  "tests/test-gls-param-gwas.R"
}
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))

Rcpp::sourceCpp(file.path(pkg_root, "src", "stage2_se_helper.cpp"))
source(file.path(pkg_root, "R", "utils.R"))
source(file.path(pkg_root, "R", "covariance.R"))
source(file.path(pkg_root, "R", "weight-functions.R"))
source(file.path(pkg_root, "R", "gwas.R"))

upper_row_major <- function(Sigma) {
  out <- numeric(nrow(Sigma) * (nrow(Sigma) + 1) / 2)
  idx <- 1L
  for (r in seq_len(nrow(Sigma))) {
    for (c in r:ncol(Sigma)) {
      out[idx] <- Sigma[r, c]
      idx <- idx + 1L
    }
  }
  out
}

manual_gls <- function(beta, Sigma, W) {
  SinvW <- solve(Sigma, W)
  SinvBeta <- solve(Sigma, beta)
  Vtheta <- solve(crossprod(W, SinvW))
  theta <- Vtheta %*% crossprod(W, SinvBeta)
  list(theta = theta, se = sqrt(diag(Vtheta)))
}

taus <- c(0.2, 0.5, 0.8)
W <- cbind(intercept = 1, slope = c(-1, 0, 1))
theta_true <- cbind(
  snp1 = c(0.10, 0.20),
  snp2 = c(-0.05, 0.15)
)
Q_slope <- W %*% theta_true + matrix(
  c(0.02, -0.01, 0.03, -0.03, 0.01, -0.02),
  nrow = 3
)
SE_tau <- matrix(0.1, nrow = 3, ncol = 2)

Sigma1 <- matrix(c(
  0.09, 0.05, 0.02,
  0.05, 0.10, 0.04,
  0.02, 0.04, 0.08
), nrow = 3, byrow = TRUE)
Sigma2 <- matrix(c(
  0.07, 0.03, 0.01,
  0.03, 0.11, 0.06,
  0.01, 0.06, 0.12
), nrow = 3, byrow = TRUE)

cov_vec <- c(upper_row_major(Sigma1), upper_row_major(Sigma2))
offsets <- c(0, length(upper_row_major(Sigma1)) * 4)

stage1 <- list(
  taus = taus,
  q_tau = rep(0, length(taus)),
  Q_slope = Q_slope,
  SE_tau = SE_tau
)

fit <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W),
  se_mode = "tau_cov",
  cov_data = list(cov_vec = cov_vec, offsets = offsets),
  return_cov = "upper"
)

expected1 <- manual_gls(Q_slope[, 1], Sigma1, W)
expected2 <- manual_gls(Q_slope[, 2], Sigma2, W)
expected_theta <- cbind(expected1$theta, expected2$theta)
expected_se <- cbind(expected1$se, expected2$se)

stopifnot(identical(fit$tau_cov_estimator, "gls"))
stopifnot(max(abs(fit$params - expected_theta)) < 1e-10)
stopifnot(max(abs(fit$SE_params - expected_se)) < 1e-10)
stopifnot(nrow(fit$param_cov_upper) == 2L)
stopifnot(ncol(fit$param_cov_upper) == 3L)

fit_ols <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W),
  se_mode = "tau_cov",
  cov_data = list(cov_vec = cov_vec, offsets = offsets),
  tau_cov_estimator = "ols"
)

stopifnot(identical(fit_ols$tau_cov_estimator, "ols"))
stopifnot(max(abs(fit_ols$params - solve(crossprod(W), t(W)) %*% Q_slope)) < 1e-10)

cat("OK: tau_cov param_gwas GLS estimates and OLS fallback match expected values\n")
