

#' Quantile-level RIF GWAS
#'
#' Runs a genome-wide association scan at user-specified quantile levels (taus)
#' using Recentered Influence Functions (RIF). This stage estimates per-tau
#' SNP slopes and their naive standard errors, independent of any parametric
#' system.
#'
#' @param Y Numeric vector (length N), phenotype.
#' @param G Numeric matrix (N x P) of SNP dosages or genotypes.
#' @param taus Numeric vector of quantile levels (default \code{seq(0.10,0.90,0.05)}).
#' @param C Optional N x K covariate matrix.
#' @param residualize_Y Logical; if TRUE, regress \code{Y} on \code{C} before quantiles
#'   (and SNPs on \code{C} too).
#' @param density_floor Positive scalar; lower bound for estimated densities
#'   \eqn{f_Y(q_tau)}.
#' @param benchmark Logical; if TRUE, include timing information.
#' @param verbose Logical; print progress messages.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{taus} – quantile levels
#'   \item \code{q_tau} – estimated baseline quantiles of Y
#'   \item \code{fhat_tau} – estimated densities at \code{q_tau}
#'   \item \code{Q_slope} – (T x P) matrix of per-tau SNP slopes
#'   \item \code{SE_tau} – (T x P) matrix of per-tau SEs
#'   \item \code{timing} – (if \code{benchmark=TRUE}) timing breakdown
#' }
#'
#' @examples
#' set.seed(1)
#' N <- 2000; P <- 20
#' taus <- seq(0.1, 0.9, 0.1)
#' G <- matrix(rbinom(N * P, 2, 0.3), N, P)
#' Y <- rnorm(N)
#' res <- quantile_gwas(Y, G, taus = taus)
#' str(res)
#'
#' @export
quantile_gwas <- function(
  Y, G,
  taus = seq(0.10, 0.90, 0.05),
  C = NULL,
  residualize_Y = FALSE,
  density_floor = 1e-8,
  benchmark = TRUE,
  verbose = TRUE
) {
  Y <- as.numeric(Y)
  G <- as.matrix(G)
  N <- length(Y); P <- ncol(G); Tt <- length(taus)
  if (nrow(G) != N) stop("nrow(G) must equal length(Y).")

  t_all0 <- proc.time()[["elapsed"]]

  # (1) Residualize Y (optional) and build RIF
  if (!is.null(C) && isTRUE(residualize_Y)) {
    if (verbose) message("Residualizing Y on C...")
    Y_star <- .residualize_on_C(Y, C)
  } else {
    Y_star <- Y
  }

  if (verbose) message("Building RIF matrix...")
  t0 <- proc.time()[["elapsed"]]
  rif <- .build_rif_matrix(Y_star, taus, density_floor = density_floor)
  RIF     <- rif$R
  q_tau   <- rif$q_tau
  fhat    <- rif$fhat_tau
  R_mean  <- colMeans(RIF)
  t1 <- proc.time()[["elapsed"]]
  time_rif <- t1 - t0

  # (2) Per-SNP tau-slopes and naive SEs
  if (verbose) message("Computing per-SNP tau-slopes...")
  t2 <- proc.time()[["elapsed"]]
  Q_slope <- matrix(NA_real_, nrow = Tt, ncol = P,
                    dimnames = list(paste0("tau", taus), colnames(G)))
  SE_tau  <- matrix(NA_real_, nrow = Tt, ncol = P,
                    dimnames = list(paste0("tau", taus), colnames(G)))

  for (j in seq_len(P)) {
    x <- G[, j]
    if (!is.null(C)) x <- .residualize_on_C(x, C)
    x_c <- x - mean(x)
    Sxx <- sum(x_c * x_c)
    if (!(is.finite(Sxx) && Sxx > 0)) next

    b <- as.numeric(crossprod(x_c, RIF) / Sxx)
    a <- R_mean - b * mean(x)
    term1 <- matrix(rep(a, each = N), nrow = N, ncol = Tt)
    term2 <- outer(x, b)
    E     <- RIF - term1 - term2

    SSE    <- colSums(E * E)
    sigma2 <- SSE / max(N - 2L, 1L)
    se_b   <- sqrt(sigma2 / Sxx)

    Q_slope[, j] <- b
    SE_tau[, j]  <- se_b
  }
  t3 <- proc.time()[["elapsed"]]
  time_gwas <- t3 - t2

  t_all1 <- proc.time()[["elapsed"]]
  timing <- list(rif_build_sec = time_rif,
                 gwas_loop_sec = time_gwas,
                 total_sec     = t_all1 - t_all0)

  out <- list(
    taus     = taus,
    q_tau    = q_tau,
    fhat_tau = fhat,
    Q_slope  = Q_slope,
    SE_tau   = SE_tau
  )
  if (benchmark) out$timing <- timing
  out
}


#' Parametric GWAS mapping
#'
#' Maps quantile-level GWAS results into effects on user-specified
#' parametric system parameters via a weight matrix \code{W}.
#'
#' @param stage1 Output list from \code{quantile_gwas()}.
#' @param transform One of:
#'   \itemize{
#'     \item \code{"custom_W"} – supply a weight matrix \code{W} in \code{transform_args}.
#'     \item \code{"two_normal"} – construct weights for a two-component normal mixture.
#'   }
#' @param transform_args Arguments for the chosen transform.
#' @param se_mode Standard error mode:
#'   \itemize{
#'     \item \code{"diagonal"} – assume independence across taus (fast).
#'     \item \code{"plugin_cor"} – plugin correlation with near-PD repair.
#'     \item \code{"dwls"} – diagonal weighted least squares (Q test calibrated).
#'   }
#' @param plugin_R Optional T x T correlation matrix for \code{"plugin_cor"}.
#'
#' @return A list with elements:
#'   \item \code{W}, \code{A}, \code{params}, \code{SE_params}, \code{Q}, \code{df}.
#'
#' @export
param_gwas <- function(
  stage1,
  transform = c("custom_W", "two_normal"),
  transform_args = list(),
  se_mode = c("diagonal", "plugin_cor", "dwls"),
  plugin_R = NULL
) {
  transform <- match.arg(transform)
  se_mode   <- match.arg(se_mode)

  taus    <- stage1$taus
  q_tau   <- stage1$q_tau
  Q_slope <- stage1$Q_slope
  SE_tau  <- stage1$SE_tau
  Tt <- length(taus)
  P  <- ncol(Q_slope)

  # (1) Construct W
  if (transform == "custom_W") {
    if (is.null(transform_args$W)) stop("Provide transform_args$W.")
    W <- as.matrix(transform_args$W)
    if (nrow(W) != Tt) stop("nrow(W) must equal length(taus).")
  } else {
    req <- c("p1","mu1","sd1","mu2","sd2","include_membership")
    missing_args <- setdiff(req, names(transform_args))
    if (length(missing_args))
      stop("Missing args for two_normal: ", paste(missing_args, collapse=", "))
    W <- make_weights_normal_mixture(
      taus = taus, q_tau = q_tau,
      p1 = transform_args$p1,
      mu1 = transform_args$mu1, sd1 = transform_args$sd1,
      mu2 = transform_args$mu2, sd2 = transform_args$sd2,
      include_membership = isTRUE(transform_args$include_membership)
    )
  }
  K <- ncol(W)
  A <- solve(crossprod(W), t(W))

  # (2) Map tau-slopes -> parameters
  params <- A %*% Q_slope
  rownames(params) <- colnames(W)
  colnames(params) <- colnames(Q_slope)

  # (3) SEs + Q statistics
  SE_params <- matrix(NA_real_, nrow = K, ncol = P,
                      dimnames = list(colnames(W), colnames(Q_slope)))
  Q_stat <- rep(NA_real_, P)

  Q_stat <- rep(NA_real_, P)

if (se_mode == "diagonal") {
  A2 <- A * A
  Vb <- SE_tau * SE_tau
  SE_params <- sqrt(A2 %*% Vb)
  
  # Q statistic (naive, assumes independence across taus)
  for (j in seq_len(P)) {
    resid_j <- Q_slope[, j] - W %*% params[, j]
    var_j   <- SE_tau[, j]^2
    Q_stat[j] <- sum((resid_j^2) / pmax(var_j, 1e-12))
  }

} else if (se_mode == "plugin_cor") {
  if (is.null(plugin_R)) {
    R_tau <- stats::cor(t(Q_slope), use = "pairwise.complete.obs")
  } else {
    R_tau <- as.matrix(plugin_R)
  }
  R_tau <- .nearPD_eig(R_tau, eps = 1e-8)

  for (j in seq_len(P)) {
    d <- SE_tau[, j]
    Dj <- diag(d, nrow = Tt, ncol = Tt)
    Sigma_j <- Dj %*% R_tau %*% Dj
    Vtheta  <- A %*% Sigma_j %*% t(A)
    Vtheta  <- .nearPD_eig(Vtheta, eps = 1e-12)
    SE_params[, j] <- sqrt(pmax(diag(Vtheta), 0))
    
    # Q statistic with full covariance
    resid_j <- Q_slope[, j] - W %*% params[, j]
    Q_stat[j] <- as.numeric(t(resid_j) %*% solve(Sigma_j, resid_j))
  }

} else if (se_mode == "dwls") {
  for (j in seq_len(P)) {
    var_j <- SE_tau[, j]^2
    Wgt   <- diag(1 / pmax(var_j, 1e-12), Tt)
    # GLS estimator (diagonal weight = DWLS)
    A_dwls <- solve(t(W) %*% Wgt %*% W, t(W) %*% Wgt)
    theta_j <- A_dwls %*% Q_slope[, j]
    params[, j] <- theta_j
    # SEs
    cov_theta <- solve(t(W) %*% Wgt %*% W)
    SE_params[, j] <- sqrt(diag(cov_theta))
    # Q statistic (DWLS)
    resid_j <- Q_slope[, j] - W %*% theta_j
    Q_stat[j] <- sum((resid_j^2) / var_j)
  }
}


  out <- list(
    W         = W,
    A         = A,
    params    = params,
    SE_params = SE_params,
    Q         = Q_stat,
    df        = Tt - K
  )
  if (se_mode == "plugin_cor") out$R_tau <- R_tau
  out
}
