#' =============================================================================
#' Generic Identifiability Diagnostics for Stage 2 Quantile GWAS
#' =============================================================================
#'
#' Given any W matrix and Sigma_tau, compute:
#' 1. How collinear are the columns (in GLS metric)?
#' 2. What's the SE inflation for each parameter due to collinearity?
#' 3. What's the best possible SE if we reparameterised optimally?
#' 4. What linear combinations are well vs poorly identified?
#'
#' =============================================================================

library(Matrix)

# Helpers for empirical Sigma_tau from Stage 1 outputs
read_stage1_taus <- function(stage1_path) {
  con <- gzfile(stage1_path, open = "rt")
  on.exit(close(con))
  header <- readLines(con, n = 1)
  if (length(header) == 0) stop("Empty stage1 file: ", stage1_path)
  cols <- strsplit(header, "\t")[[1]]
  tau_cols <- grep("^beta_tau_", cols, value = TRUE)
  if (length(tau_cols) == 0) tau_cols <- grep("^se_tau_", cols, value = TRUE)
  if (length(tau_cols) == 0) stop("No beta_tau_*/se_tau_* columns found in: ", stage1_path)
  taus <- as.numeric(sub("^(beta|se)_tau_", "", tau_cols))
  if (any(is.na(taus))) stop("Failed to parse tau values from header in: ", stage1_path)
  sort(unique(taus))
}

upper_to_mat <- function(upper, T) {
  mat <- matrix(0, T, T)
  mat[upper.tri(mat, diag = TRUE)] <- upper
  mat <- mat + t(mat) - diag(diag(mat))
  mat
}

mean_cov_from_cov_gz <- function(cov_path, T, max_snps = 20000L) {
  n_elem <- as.integer(T * (T + 1) / 2)
  scale_path <- sub("\\.cov\\.gz$", ".cov.scale.gz", cov_path)
  con <- gzfile(cov_path, open = "rb")
  on.exit(close(con))
  con_scale <- NULL
  if (file.exists(scale_path)) {
    con_scale <- gzfile(scale_path, open = "rb")
    on.exit(close(con_scale), add = TRUE)
  }

  sum_upper <- numeric(n_elem)
  n_read <- 0L
  repeat {
    if (is.null(con_scale)) {
      vec <- readBin(con, what = "numeric", n = n_elem, size = 4, endian = "little")
      if (length(vec) != n_elem) break
    } else {
      vec_q <- readBin(con, what = "integer", n = n_elem, size = 1, signed = TRUE, endian = "little")
      if (length(vec_q) != n_elem) break
      sc <- readBin(con_scale, what = "numeric", n = 1L, size = 4, endian = "little")
      if (length(sc) != 1L) stop("Scale file ended early: ", scale_path)
      vec <- as.numeric(vec_q) * sc
    }
    sum_upper <- sum_upper + vec
    n_read <- n_read + 1L
    if (n_read >= max_snps) break
  }

  if (n_read == 0L) stop("No covariance entries read from: ", cov_path)
  list(mean_upper = sum_upper / n_read, n_read = n_read)
}

Sigma_tau_from_stage1_cov <- function(stage1_path, cov_path, max_snps = 20000L) {
  taus <- read_stage1_taus(stage1_path)
  T <- length(taus)
  cov_mean <- mean_cov_from_cov_gz(cov_path, T, max_snps = max_snps)
  Sigma_tau <- upper_to_mat(cov_mean$mean_upper, T)
  # Stabilize to be positive definite if needed (empirical cov can be noisy)
  if (min(eigen(Sigma_tau, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    Sigma_tau <- as.matrix(nearPD(Sigma_tau, corr = FALSE)$mat)
  }
  list(Sigma_tau = Sigma_tau, taus = taus, n_read = cov_mean$n_read)
}

#' Comprehensive identifiability diagnostics
#'
#' @param W Weight matrix (K quantiles × P parameters)
#' @param Sigma_tau Covariance matrix of tau-slopes (K × K)
#' @param param_names Optional names for parameters
#'
#' @return List with diagnostic information
#'
#' @examples
#' taus <- seq(0.1, 0.9, by = 0.1)
#' W <- cbind(additive = 1, multiplicative = 25 + 4 * qnorm(taus))
#' Sigma_tau <- make_Sigma_tau(taus, n = 10000)
#' diag <- diagnose_identifiability(W, Sigma_tau)
 diagnose_identifiability <- function(W, Sigma_tau, param_names = NULL) {

  K <- nrow(W)  # Number of quantiles
  P <- ncol(W)  # Number of parameters

  if (is.null(param_names)) {
    param_names <- colnames(W)
    if (is.null(param_names)) {
      param_names <- paste0("theta_", seq_len(P))
    }
  }

  # ===========================================================================
  # 1. Basic GLS quantities
  # ===========================================================================

  Sigma_inv <- solve(Sigma_tau)

  # Information matrix (precision of theta estimates)
  information <- t(W) %*% Sigma_inv %*% W

  # Covariance of theta estimates
  covariance <- solve(information)

  # Standard errors
  se <- sqrt(diag(covariance))
  names(se) <- param_names

  # ===========================================================================
  # 2. Eigendecomposition - what directions are identifiable?
  # ===========================================================================

  eig <- eigen(information, symmetric = TRUE)
  eigenvalues <- eig$values
  eigenvectors <- eig$vectors
  rownames(eigenvectors) <- param_names
  colnames(eigenvectors) <- paste0("PC", seq_len(P))

  # Condition number
  condition_number <- max(eigenvalues) / min(eigenvalues)

  # SE in each eigendirection
  se_eigen <- 1 / sqrt(eigenvalues)
  names(se_eigen) <- paste0("PC", seq_len(P))

  # ===========================================================================
  # 3. SE inflation factors - how much worse is each parameter due to collinearity?
  # ===========================================================================

  # "Oracle" SE: what we'd get if we only estimated this parameter
  # (i.e., the other parameters were known/fixed)
  se_oracle <- numeric(P)
  for (j in seq_len(P)) {
    # Information for parameter j alone
    info_j <- t(W[, j]) %*% Sigma_inv %*% W[, j]
    se_oracle[j] <- 1 / sqrt(info_j)
  }
  names(se_oracle) <- param_names

  # SE inflation factor
  se_inflation <- se / se_oracle
  names(se_inflation) <- param_names

  # Variance inflation factor (VIF) - traditional measure
  # VIF_j = Var(theta_hat_j) / Var(theta_hat_j | others fixed)
  #       = se_inflation^2
  vif <- se_inflation^2
  names(vif) <- param_names

  # ===========================================================================
  # 4. GLS correlation/cosine structure
  # ===========================================================================

  # Correlation between parameter estimates
  D <- diag(1/se)
  cor_estimates <- D %*% covariance %*% D
  rownames(cor_estimates) <- colnames(cor_estimates) <- param_names

  # GLS cosine similarities between W columns
  # (This measures collinearity in the Sigma^{-1} metric)
  gls_norms <- sqrt(diag(information))
  gls_cosine <- information / outer(gls_norms, gls_norms)
  rownames(gls_cosine) <- colnames(gls_cosine) <- param_names

  # ===========================================================================
  # 5. Best possible SE after optimal reparameterisation
  # ===========================================================================

  # The best SE for detecting *any* effect is along PC1
  best_se <- min(se_eigen)  # = 1/sqrt(max eigenvalue)

  # The worst SE (for the hardest-to-identify direction) is along PC_last
  worst_se <- max(se_eigen)  # = 1/sqrt(min eigenvalue)

  # ===========================================================================
  # 6. Optimal orthogonal reparameterisation
  # ===========================================================================

  # W_orthogonal: columns are GLS-orthogonal
  # Use Cholesky of Sigma_inv for whitening
  Sigma_inv_sqrt <- chol(Sigma_inv)  # Upper triangular
  W_whitened <- Sigma_inv_sqrt %*% W

  # QR decomposition of whitened W
  qr_result <- qr(W_whitened)
  Q <- qr.Q(qr_result)  # Orthonormal in standard metric
  R <- qr.R(qr_result)  # Transformation matrix

  # Orthogonalised W (back in original space)
  W_orthogonal <- solve(Sigma_inv_sqrt) %*% Q

  # Transformation: theta_original = R %*% theta_orthogonal
  # (R is upper triangular, so theta_orth[1] is a combo of all original params, etc.)

  # SE in orthogonal parameterisation
  info_orth <- t(W_orthogonal) %*% Sigma_inv %*% W_orthogonal
  se_orthogonal <- 1 / sqrt(diag(info_orth))

  # ===========================================================================
  # 7. Summary statistics
  # ===========================================================================

  # Maximum possible SE reduction (if we could perfectly orthogonalise
  # and test along the best direction)
  max_se_reduction <- max(se) / best_se

  # Geometric mean SE inflation
  mean_se_inflation <- exp(mean(log(se_inflation)))

  # ===========================================================================
  # Return everything
  # ===========================================================================

  result <- list(
    # Basic quantities
    information = information,
    covariance = covariance,
    se = se,

    # Eigenstructure
    eigenvalues = eigenvalues,
    eigenvectors = eigenvectors,
    condition_number = condition_number,
    se_eigen = se_eigen,

    # Inflation factors
    se_oracle = se_oracle,
    se_inflation = se_inflation,
    vif = vif,

    # Correlation structure
    cor_estimates = cor_estimates,
    gls_cosine = gls_cosine,

    # Optimal reparameterisation
    W_orthogonal = W_orthogonal,
    R_transform = R,  # theta_original = R %*% theta_orthogonal
    se_orthogonal = se_orthogonal,

    # Summary
    best_se = best_se,
    worst_se = worst_se,
    max_se_reduction = max_se_reduction,
    mean_se_inflation = mean_se_inflation,

    # Metadata
    param_names = param_names,
    n_quantiles = K,
    n_params = P
  )

  class(result) <- "identifiability_diagnostic"
  return(result)
}


#' Print method for identifiability diagnostics
print.identifiability_diagnostic <- function(x, ...) {

  cat("=======================================================================\n")
  cat("IDENTIFIABILITY DIAGNOSTICS\n")
  cat("=======================================================================\n\n")

  cat(sprintf("Parameters: %d    Quantiles: %d\n\n", x$n_params, x$n_quantiles))

  # Condition number interpretation
  cat("CONDITION NUMBER: ", sprintf("%.1f", x$condition_number), "\n")
  if (x$condition_number < 10) {
    cat("  → Well-conditioned. Parameters are clearly distinguishable.\n\n")
  } else if (x$condition_number < 100) {
    cat("  → Moderately ill-conditioned. Some SE inflation expected.\n\n")
  } else if (x$condition_number < 1000) {
    cat("  → Ill-conditioned. Substantial SE inflation.\n\n")
  } else {
    cat("  → Severely ill-conditioned. Consider reparameterisation.\n\n")
  }

  # SE table
  cat("STANDARD ERRORS:\n")
  cat(sprintf("  %-20s %12s %12s %12s\n", "Parameter", "Current SE", "Oracle SE", "Inflation"))
  cat("  ", strrep("-", 58), "\n", sep = "")
  for (j in seq_along(x$param_names)) {
    cat(sprintf("  %-20s %12.5f %12.5f %12.2fx\n",
                x$param_names[j], x$se[j], x$se_oracle[j], x$se_inflation[j]))
  }
  cat("\n")

  # Correlation matrix
  cat("CORRELATION OF ESTIMATES:\n")
  print(round(x$cor_estimates, 3))
  cat("\n")

  # Eigenvectors
  cat("IDENTIFIABLE DIRECTIONS (eigenvectors of information matrix):\n")
  for (k in seq_len(x$n_params)) {
    direction <- paste(sprintf("%.3f×%s", x$eigenvectors[, k], x$param_names),
                       collapse = " + ")
    cat(sprintf("  PC%d (SE = %.5f): %s\n", k, x$se_eigen[k], direction))
  }
  cat("\n")

  # Summary
  cat("SUMMARY:\n")
  cat(sprintf("  Best achievable SE (along PC1):     %.5f\n", x$best_se))
  cat(sprintf("  Worst SE (along PC%d):              %.5f\n", x$n_params, x$worst_se))
  cat(sprintf("  Maximum potential SE reduction:     %.1fx\n", x$max_se_reduction))
  cat(sprintf("  Mean SE inflation factor:           %.2fx\n", x$mean_se_inflation))

  invisible(x)
}


#' Suggest optimal reparameterisation
#'
#' Given diagnostics, suggest how to reparameterise for better identifiability
suggest_reparameterisation <- function(diag) {

  cat("=======================================================================\n")
  cat("REPARAMETERISATION SUGGESTIONS\n")
  cat("=======================================================================\n\n")

  if (diag$condition_number < 10) {
    cat("Condition number is low. No reparameterisation needed.\n\n")
    return(invisible(NULL))
  }

  # Find which parameters are most problematic
  worst_param <- which.max(diag$se_inflation)

  cat(sprintf("Most inflated parameter: %s (%.1fx inflation)\n\n",
              diag$param_names[worst_param], diag$se_inflation[worst_param]))

  # Look at GLS cosine to find what it's collinear with
  cosines <- diag$gls_cosine[worst_param, ]
  cosines[worst_param] <- 0  # Exclude self
  most_collinear <- which.max(abs(cosines))

  cat(sprintf("Most collinear with: %s (GLS cosine = %.3f)\n\n",
              diag$param_names[most_collinear], cosines[most_collinear]))

  cat("SUGGESTED ACTIONS:\n\n")

  cat("1. ORTHOGONALISE: Replace W with W_orthogonal from diagnostics.\n")
  cat("   New parameters are linear combinations of original:\n")
  cat("   theta_original = R %*% theta_orthogonal\n")
  cat("   where R is:\n")
  print(round(diag$R_transform, 3))
  cat("\n")

  cat("2. CENTRE: If one column of W has a large constant component,\n")
  cat("   subtract the mean to remove correlation with the intercept.\n\n")

  cat("3. REPORT IDENTIFIABLE QUANTITIES: Instead of original parameters,\n")
  cat("   report the well-identified eigendirections:\n")
  for (k in seq_len(min(2, diag$n_params))) {
    direction <- paste(sprintf("%.3f×%s", diag$eigenvectors[, k], diag$param_names),
                       collapse = " + ")
    cat(sprintf("   PC%d = %s (SE = %.5f)\n", k, direction, diag$se_eigen[k]))
  }
  cat("\n")

  invisible(diag)
}


#' Compute theoretical Sigma_tau for quantile regression
#'
#' @param taus Vector of quantile levels
#' @param n Sample size
#' @param density_at_quantile Function returning density at each quantile
#'        Default assumes standard normal.
make_Sigma_tau <- function(taus, n, density_at_quantile = NULL) {

  K <- length(taus)

  if (is.null(density_at_quantile)) {
    density_at_quantile <- function(tau) dnorm(qnorm(tau))
  }

  Sigma <- matrix(NA, K, K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      f_i <- density_at_quantile(taus[i])
      f_j <- density_at_quantile(taus[j])
      Sigma[i, j] <- (min(taus[i], taus[j]) - taus[i] * taus[j]) / (n * f_i * f_j)
    }
  }

  return(Sigma)
}


# =============================================================================
# EXAMPLES
# =============================================================================

if (isTRUE(getOption("run_identifiability_examples", TRUE))) {  # Run examples

  cat("\n\n")
  cat("#######################################################################\n")
  cat("# EXAMPLE 1: Single Normal - Raw vs Centred Parameterisation\n")
  cat("#######################################################################\n\n")

  taus <- seq(0.1, 0.9, by = 0.1)
  mu <- 25
  sigma <- 4
  n <- 10000

  Sigma_tau <- make_Sigma_tau(taus, n)

  # Raw parameterisation (problematic)
  W_raw <- cbind(
    additive = rep(1, length(taus)),
    multiplicative = mu + sigma * qnorm(taus)
  )

  cat("--- RAW PARAMETERISATION (multiplicative = Q(tau)) ---\n\n")
  diag_raw <- diagnose_identifiability(W_raw, Sigma_tau)
  print(diag_raw)

  # Centred parameterisation (fixed)
  W_centred <- cbind(
    additive = rep(1, length(taus)),
    multiplicative = sigma * qnorm(taus)
  )

  cat("\n--- CENTRED PARAMETERISATION (multiplicative = Q(tau) - mu) ---\n\n")
  diag_centred <- diagnose_identifiability(W_centred, Sigma_tau)
  print(diag_centred)

  # Summary comparison
  cat("\n--- COMPARISON ---\n")
  cat(sprintf("Condition number:  Raw = %.1f,  Centred = %.1f\n",
              diag_raw$condition_number, diag_centred$condition_number))
  cat(sprintf("SE(additive):      Raw = %.5f,  Centred = %.5f\n",
              diag_raw$se[1], diag_centred$se[1]))
  cat(sprintf("SE improvement:    %.1fx\n",
              diag_raw$se[1] / diag_centred$se[1]))


  cat("\n\n")
  cat("#######################################################################\n")
  cat("# EXAMPLE 2: Two-Component Mixture (well-separated)\n")
  cat("#######################################################################\n\n")

  # Numerical derivatives for mixture W matrix
  mixture_quantile <- function(tau, p1, mu1, sigma1, mu2, sigma2) {
    objective <- function(x) {
      p1 * pnorm(x, mu1, sigma1) + (1 - p1) * pnorm(x, mu2, sigma2) - tau
    }
    uniroot(objective, c(-100, 200))$root
  }

  p1 <- 0.6
  mu1 <- 22
  sigma1 <- 2
  mu2 <- 28
  sigma2 <- 4
  eps <- 1e-5

  W_mixture <- matrix(NA, length(taus), 3)
  colnames(W_mixture) <- c("p1", "mu1", "mu2")

  base_q <- sapply(taus, mixture_quantile, p1, mu1, sigma1, mu2, sigma2)

  for (i in seq_along(taus)) {
    # dQ/dp1
    q_p <- mixture_quantile(taus[i], p1 + eps, mu1, sigma1, mu2, sigma2)
    W_mixture[i, 1] <- (q_p - base_q[i]) / eps

    # dQ/dmu1
    q_m1 <- mixture_quantile(taus[i], p1, mu1 + eps, sigma1, mu2, sigma2)
    W_mixture[i, 2] <- (q_m1 - base_q[i]) / eps

    # dQ/dmu2
    q_m2 <- mixture_quantile(taus[i], p1, mu1, sigma1, mu2 + eps, sigma2)
    W_mixture[i, 3] <- (q_m2 - base_q[i]) / eps
  }

  cat("--- MIXTURE MODEL (p1, mu1, mu2) ---\n\n")
  diag_mixture <- diagnose_identifiability(W_mixture, Sigma_tau)
  print(diag_mixture)

  cat("\n")
  suggest_reparameterisation(diag_mixture)

  cat("\n\n")
  cat("#######################################################################\n")
  cat("# EXAMPLE 3: Two-Component Mixture (overlapped)\n")
  cat("#######################################################################\n\n")

  p1b <- 0.5
  mu1b <- 24
  sigma1b <- 3
  mu2b <- 26
  sigma2b <- 3

  W_mixture_b <- matrix(NA, length(taus), 3)
  colnames(W_mixture_b) <- c("p1", "mu1", "mu2")

  base_q_b <- sapply(taus, mixture_quantile, p1b, mu1b, sigma1b, mu2b, sigma2b)

  for (i in seq_along(taus)) {
    q_p <- mixture_quantile(taus[i], p1b + eps, mu1b, sigma1b, mu2b, sigma2b)
    W_mixture_b[i, 1] <- (q_p - base_q_b[i]) / eps

    q_m1 <- mixture_quantile(taus[i], p1b, mu1b + eps, sigma1b, mu2b, sigma2b)
    W_mixture_b[i, 2] <- (q_m1 - base_q_b[i]) / eps

    q_m2 <- mixture_quantile(taus[i], p1b, mu1b, sigma1b, mu2b + eps, sigma2b)
    W_mixture_b[i, 3] <- (q_m2 - base_q_b[i]) / eps
  }

  cat("--- MIXTURE MODEL (p1, mu1, mu2) OVERLAPPED ---\n\n")
  diag_mixture_b <- diagnose_identifiability(W_mixture_b, Sigma_tau)
  print(diag_mixture_b)


  cat("\n\n")
  cat("#######################################################################\n")
  cat("# EXAMPLE 3: vQTL (mean + variance effect)\n")
  cat("#######################################################################\n\n")

  W_vqtl <- cbind(
    beta_mu = rep(1, length(taus)),
    beta_sigma2 = sigma * qnorm(taus)
  )
  diag_vqtl <- diagnose_identifiability(W_vqtl, Sigma_tau)
  print(diag_vqtl)

  cat("\n\n")
  cat("#######################################################################\n")
  cat("# EXAMPLE 4: Effect of mu/sigma ratio\n")
  cat("#######################################################################\n\n")

  cat(sprintf("%-12s %15s %15s %15s\n",
              "mu/sigma", "Condition #", "SE inflation", "Correlation"))
  cat(strrep("-", 60), "\n")

  for (ratio in c(0, 1, 2, 5, 10, 25)) {
    mu_test <- ratio * sigma
    W_test <- cbind(
      additive = rep(1, length(taus)),
      multiplicative = mu_test + sigma * qnorm(taus)
    )
    diag_test <- diagnose_identifiability(W_test, Sigma_tau)

    cat(sprintf("%-12.0f %15.1f %15.2fx %15.3f\n",
                ratio,
                diag_test$condition_number,
                diag_test$se_inflation[1],
                diag_test$cor_estimates[1, 2]))
  }

  cat("\nNote: At mu/sigma = 0, the parameterisations are equivalent.\n")
  cat("As mu/sigma increases, the raw parameterisation becomes increasingly\n")
  cat("ill-conditioned because Q(tau) ≈ mu × constant + small deviation.\n")
}

if (FALSE) {
  # 1. Construct your W matrix (whatever parameterisation you're using)
  W <- make_weights_whatever(taus, params)

  # 2. Get Sigma_tau from stage 1 (or theoretical)
  Sigma_tau <- estimate_Sigma_tau(stage1_results)

  # 3. Run diagnostic
  diag <- diagnose_identifiability(W, Sigma_tau)

  # 4. Check SE inflation - this is your "price"
  print(diag$se_inflation)
  # If any parameter has inflation > 2, consider reparameterising

  # 5. Check if centring helps
  W_centred <- centre_W_columns(W)  # Your choice of how to centre
  diag_centred <- diagnose_identifiability(W_centred, Sigma_tau)
  print(diag_centred$se_inflation)
}
