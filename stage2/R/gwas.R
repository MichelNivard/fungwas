#' Quantile-level RIF GWAS
#'
#' Runs a genome-wide association scan at user-specified quantile levels (taus)
#' using Recentered Influence Functions (RIF). Supports two covariate modes:
#' - "residualize": partial out C from SNPs (but never from Y).
#' - "include": include C directly in the regression alongside each SNP.
#'
#' @param Y Numeric vector (length N), phenotype.
#' @param G Numeric matrix (N x P) of SNP dosages or genotypes.
#' @param taus Numeric vector of quantile levels.
#' @param C Optional N x K covariate matrix.
#' @param covar_mode Either "residualize" (default) or "include".
#'   - "residualize": partial out C from G before regression.
#'   - "include": include C directly in regression.
#' @param density_floor Positive scalar; lower bound for estimated densities.
#' @param benchmark Logical; if TRUE, include timing info.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with taus, q_tau, fhat_tau, Q_slope, SE_tau, timing.
#' 
#' @export
quantile_gwas <- function(
  Y, G,
  taus = seq(0.10, 0.90, 0.05),
  C = NULL,
  covar_mode = c("residualize", "include"),
  density_floor = 1e-8,
  benchmark = TRUE,
  verbose = TRUE
) {
  covar_mode <- match.arg(covar_mode)
  Y <- as.numeric(Y)
  G <- as.matrix(G)
  N <- length(Y); P <- ncol(G); Tt <- length(taus)
  if (nrow(G) != N) stop("nrow(G) must equal length(Y).")

  t_all0 <- proc.time()[["elapsed"]]

  # (1) Build RIF on raw Y
  if (verbose) message("Building RIF matrix on raw Y...")
  t0 <- proc.time()[["elapsed"]]
  rif <- .build_rif_matrix(Y, taus, density_floor = density_floor)
  RIF     <- rif$R
  q_tau   <- rif$q_tau
  fhat    <- rif$fhat_tau
  R_mean  <- colMeans(RIF)
  t1 <- proc.time()[["elapsed"]]
  time_rif <- t1 - t0

  # (2) SNP regressions
  if (verbose) message("Computing per-SNP tau-slopes...")
  t2 <- proc.time()[["elapsed"]]
  Q_slope <- matrix(NA_real_, nrow = Tt, ncol = P,
                    dimnames = list(paste0("tau", taus), colnames(G)))
  SE_tau  <- matrix(NA_real_, nrow = Tt, ncol = P,
                    dimnames = list(paste0("tau", taus), colnames(G)))

  if (covar_mode == "residualize") {
    # Partial out C from G (never from Y)
    if (!is.null(C)) {
      if (verbose) message("Residualizing SNPs on covariates...")
      for (j in seq_len(P)) {
        G[, j] <- .residualize_on_C(G[, j], C)
      }
    }
    # Univariate regression of RIF on each SNP
    for (j in seq_len(P)) {
      x <- G[, j]
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

  } else if (covar_mode == "include") {
    print("Currently covar_mode = 'include' is still in development")
    break

    # Joint regression RIF ~ SNP + C
    if (is.null(C)) stop("Must provide C if covar_mode = 'include'.")
    C <- as.matrix(C)

    for (j in seq_len(P)) {
      X <- cbind(G[, j], C)  # SNP + covariates
      fit <- lm.fit(X, RIF)  # multi-response regression (N x T outcome)
      b <- fit$coefficients[1, ]      # SNP slope
      Q_slope[, j] <- b

      # compute SEs for SNP slope (per tau)
      resid <- RIF - X %*% fit$coefficients
      sigma2 <- colSums(resid^2) / (N - ncol(X))
      XtXinv <- chol2inv(chol(crossprod(X)))
      se_b <- sqrt(sigma2 * XtXinv[1, 1])
      SE_tau[, j] <- se_b
    }
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
#'     \item \code{"tau_cov"} – use per-SNP jackknife covariance from Stage 1 (most accurate).
#'   }
#' @param plugin_R Optional T x T correlation matrix for \code{"plugin_cor"}.
#' @param cov_data For \code{se_mode = "tau_cov"}: either a list with \code{cov_vec} 
#'   (numeric vector of covariance data) and \code{offsets} (byte offsets per SNP),
#'   or NULL to use diagonal fallback.
#'
#' @return A list with the following elements:
#' \item{\code{W}}{Weight matrix used in the parametric mapping.}
#' \item{\code{A}}{Mapping matrix (tau-by-parameter).}
#' \item{\code{params}}{Estimated parameter values.}
#' \item{\code{SE_params}}{Delta-method standard errors of parameter estimates.}
#' \item{\code{Q}}{Vector of Q statistics per SNP (if computed).}
#' \item{\code{df}}{Degrees of freedom for Q statistic.}
#' 
#' @export
param_gwas <- function(
  stage1,
  transform = c("custom_W", "two_normal"),
  transform_args = list(),
  se_mode = c("diagonal", "plugin_cor", "dwls", "tau_cov"),
  plugin_R = NULL,
  cov_data = NULL
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

} else if (se_mode == "tau_cov") {
  # Use per-SNP jackknife covariance from Stage 1
  if (is.null(cov_data)) {
    warning("tau_cov mode requires cov_data; falling back to diagonal")
    A2 <- A * A
    Vb <- SE_tau * SE_tau
    SE_params <- sqrt(A2 %*% Vb)
  } else {
    # Use C++ helper for fast computation
    se_mat <- compute_calibrated_se(
      cov_vec = cov_data$cov_vec,
      offsets = cov_data$offsets,
      A = A,
      K = Tt
    )
    SE_params <- t(se_mat)
    rownames(SE_params) <- colnames(W)
  }
  
  # Q statistic (using diagonal approximation)
  for (j in seq_len(P)) {
    resid_j <- Q_slope[, j] - W %*% params[, j]
    var_j   <- SE_tau[, j]^2
    Q_stat[j] <- sum((resid_j^2) / pmax(var_j, 1e-12))
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


#' Parametric GWAS from Stage 1 files
#'
#' Convenience wrapper for \code{param_gwas()} that reads Stage 1 output files
#' directly, including covariance data for \code{se_mode = "tau_cov"}.
#'
#' @param stage1_file Path to Stage 1 TSV file (.stage1.tsv.gz)
#' @param W Weight matrix (T x K)
#' @param se_mode SE computation mode (see \code{param_gwas})
#' @param cov_file Path to covariance file (.cov.gz) for tau_cov mode
#' @param rtau_file Path to R_tau correlation matrix (.Rtau.tsv) for plugin_cor mode
#'
#' @return A data.table with SNP info and parameter estimates
#' 
#' @export
param_gwas_from_file <- function(
  stage1_file,
  W,
  se_mode = c("diagonal", "plugin_cor", "dwls", "tau_cov"),
  cov_file = NULL,
  rtau_file = NULL
) {
  se_mode <- match.arg(se_mode)
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required for param_gwas_from_file")
  }
  
  # Read Stage 1
  message("Reading Stage 1 file: ", stage1_file)
  stage1_dt <- data.table::fread(stage1_file)
  
  # Extract beta/SE columns
  beta_cols <- grep("^beta_tau_", names(stage1_dt), value = TRUE)
  se_cols <- grep("^se_tau_", names(stage1_dt), value = TRUE)
  
  taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
  
  Q_slope <- t(as.matrix(stage1_dt[, ..beta_cols]))
  SE_tau <- t(as.matrix(stage1_dt[, ..se_cols]))
  
  # Build stage1 list
  stage1 <- list(
    taus = taus,
    q_tau = rep(0, length(taus)),  # Not used for param_gwas
    Q_slope = Q_slope,
    SE_tau = SE_tau
  )
  
  # Load covariance if needed
  cov_data <- NULL
  if (se_mode == "tau_cov" && !is.null(cov_file)) {
    message("Loading covariance file: ", cov_file)
    
    # Read binary covariance
    con <- gzfile(cov_file, "rb")
    cov_vec <- readBin(con, "numeric", n = 2e9, size = 4)
    close(con)
    
    # Offsets: assume sequential storage
    # T*(T+1)/2 floats per SNP, 4 bytes each
    n_cov_elements <- length(taus) * (length(taus) + 1) / 2
    n_snps <- nrow(stage1_dt)
    offsets <- seq(0, by = n_cov_elements * 4, length.out = n_snps)
    
    cov_data <- list(cov_vec = cov_vec, offsets = offsets)
  }
  
  # Load R_tau if needed
  plugin_R <- NULL
  if (se_mode == "plugin_cor" && !is.null(rtau_file)) {
    message("Loading R_tau file: ", rtau_file)
    rtau_dt <- data.table::fread(rtau_file)
    plugin_R <- as.matrix(rtau_dt)
  }
  
  # Run param_gwas
  result <- param_gwas(
    stage1 = stage1,
    transform = "custom_W",
    transform_args = list(W = W),
    se_mode = se_mode,
    plugin_R = plugin_R,
    cov_data = cov_data
  )
  
  # Build output data.table
  out_dt <- stage1_dt[, .(snp_id, chr, bp, a1, a2)]
  
  # Add parameters and SEs
  for (i in seq_len(nrow(result$params))) {
    param_name <- rownames(result$params)[i]
    out_dt[[param_name]] <- result$params[i, ]
    out_dt[[paste0(param_name, "_se")]] <- result$SE_params[i, ]
  }
  
  out_dt$Q <- result$Q
  out_dt$Q_pval <- stats::pchisq(result$Q, df = result$df, lower.tail = FALSE)
  
  out_dt
}

