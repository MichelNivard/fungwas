# ============================================================
# Test: fit_multi_models_tau_cov basic correctness
# ============================================================

suppressPackageStartupMessages(library(Rcpp))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- sub("^--file=", "", file_arg)
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))

Rcpp::sourceCpp(file.path(pkg_root, "src", "models.cpp"))

# Two SNPs, one tau
beta_stage1 <- matrix(c(0.5, -1.0), nrow = 2, ncol = 1)

# Covariance: 1x1 per SNP
cov_vec <- c(0.25, 0.50)  # variances for each SNP
# Offsets are in bytes (float32 = 4 bytes)
offsets <- c(0, 4)

W_list <- list(add = matrix(1, nrow = 1, ncol = 1))
N <- 100

res <- fit_multi_models_tau_cov(beta_stage1, cov_vec, offsets, W_list, N, n_threads = 1)

# With W = 1, theta = beta, resid = 0, so Q = 0
expected_aic <- 2
expected_bic <- log(N)

stopifnot(all(abs(res[, 1] - expected_aic) < 1e-8))
stopifnot(all(abs(res[, 2] - expected_bic) < 1e-8))

cat("OK: tau_cov multi-model AIC/BIC matches expected values\n")
