
.kde_f_at_q <- function(y, q) {
  bw <- suppressWarnings(stats::bw.nrd0(y))
  if (!is.finite(bw) || bw < 1e-8) bw <- max(stats::sd(y), 1e-3)
  mean(dnorm((y - q) / bw)) / bw
}

.build_rif_matrix <- function(y, taus, density_floor = 1e-8) {
  q    <- as.numeric(quantile(y, probs = taus, type = 8))
  fhat <- vapply(q, function(qt) max(.kde_f_at_q(y, qt), density_floor), numeric(1))
  N <- length(y); Tt <- length(taus)
  R <- matrix(NA_real_, nrow = N, ncol = Tt)
  for (t in seq_along(taus)) {
    ind <- as.numeric(y <= q[t])
    R[, t] <- q[t] + (taus[t] - ind) / fhat[t]
  }
  colnames(R) <- paste0("tau", taus)
  list(R = R, q_tau = q, fhat_tau = fhat)
}

.residualize_on_C <- function(v, C) {
  # Fast OLS residuals of v ~ C (adds intercept internally via centering trick)
  if (is.null(C)) return(v)
  C <- as.matrix(C)
  # Add intercept column
  X <- cbind(Intercept = 1, C)
  qrx <- qr(X)
  # fitted = X %*% (X^+ v)
  coef <- qr.coef(qrx, v)
  v - as.numeric(X %*% coef)
}


.nearPD_eig <- function(S, eps = 1e-8) {
  # Symmetrize, eigen clip, small jitter if needed
  S <- 0.5 * (S + t(S))
  ev <- eigen(S, symmetric = TRUE)
  vals <- pmax(ev$values, eps)
  S_pd <- ev$vectors %*% (vals * t(ev$vectors))
  0.5 * (S_pd + t(S_pd))
}

