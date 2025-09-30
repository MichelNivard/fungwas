#' Build tau-by-K weights for a 2-component Normal mixture
#'
#' Constructs a weight matrix \eqn{W} that maps RIF tau-slopes \eqn{b(\tau)}
#' into component-parameter effects under a two-Normal mixture baseline.
#'
#' Let \eqn{Y \sim p_1 N(\mu_1, \sigma_1^2) + (1-p_1) N(\mu_2, \sigma_2^2)}.
#' For baseline unconditional quantiles \eqn{q_\tau} and mixture pdf
#' \eqn{f(q_\tau)}, the default columns of \eqn{W} correspond to component means:
#' \deqn{W_1(\tau) = \frac{p_1 f_1(q_\tau)}{f(q_\tau)}, \quad
#'       W_2(\tau) = \frac{(1-p_1) f_2(q_\tau)}{f(q_\tau)}}
#'
#' If \code{include_membership = TRUE}, a first column is added for membership
#' (log-odds) perturbation:
#' \deqn{W_\gamma(\tau) = - \frac{p_1 (1-p_1)}{f(q_\tau)}
#'       \left\{F_1(q_\tau) - F_2(q_\tau)\right\}}
#'
#' @param taus Numeric vector of quantile levels.
#' @param q_tau Numeric vector of baseline quantiles at \code{taus} (\code{type=8} recommended).
#' @param p1 Proportion of class 1 in (0,1); class 2 is \code{1 - p1}.
#' @param mu1,sd1 Mean and standard deviation of component 1.
#' @param mu2,sd2 Mean and standard deviation of component 2.
#' @param include_membership Logical; include a first column \code{"gamma"} if TRUE.
#' @param tiny Small positive floor for stabilizing divisions (default: mixture pdf clamped to \code{tiny}).
#'
#' @return A numeric matrix (\code{T x K}) with column names:
#'   \itemize{
#'     \item If \code{include_membership = FALSE}: \code{c("beta_1","beta_2")}
#'     \item Else: \code{c("gamma","beta_1","beta_2")}
#'   }
#'
#' @examples
#' taus <- seq(0.10, 0.90, by = 0.05)
#' y <- rnorm(2000, 2, 1)
#' q_tau <- as.numeric(quantile(y, taus, type = 8))
#' W <- make_weights_normal_mixture(
#'   taus, q_tau,
#'   p1 = 0.5, mu1 = 1.2, sd1 = 0.45,
#'   mu2 = 3.0, sd2 = 0.7,
#'   include_membership = TRUE
#' )
#' dim(W); colnames(W)
#'
#' @export
make_weights_normal_mixture <- function(
  taus,
  q_tau,
  p1,
  mu1, sd1,
  mu2, sd2,
  include_membership = FALSE,
  tiny = 1e-12
) {
  stopifnot(length(taus) == length(q_tau))
  Tt <- length(taus)

  F1 <- pnorm(q_tau, mean = mu1, sd = sd1)
  f1 <- dnorm(q_tau, mean = mu1, sd = sd1)
  F2 <- pnorm(q_tau, mean = mu2, sd = sd2)
  f2 <- dnorm(q_tau, mean = mu2, sd = sd2)

  f  <- pmax(p1 * f1 + (1 - p1) * f2, tiny)

  W1 <- p1 * f1 / f
  W2 <- (1 - p1) * f2 / f

  if (isTRUE(include_membership)) {
    Wp <- -(p1 * (1 - p1) / f) * (F1 - F2)
    W  <- cbind(Wp, W1, W2)
    colnames(W) <- c("gamma", "beta_1", "beta_2")
  } else {
    W  <- cbind(W1, W2)
    colnames(W) <- c("beta_1", "beta_2")
  }
  rownames(W) <- paste0("tau", taus)
  storage.mode(W) <- "double"
  W
}


#' Build tau-by-2 weights for Normal mean/variance perturbations
#'
#' Constructs a weight matrix W that maps RIF tau-slopes b(tau) into
#' effects on the mean and variance of a Normal baseline phenotype.
#'
#' Baseline: Y ~ N(mu, sigma^2).
#' For tau-quantiles q_tau = mu + sigma * z_tau, z_tau = Phi^{-1}(tau):
#' \deqn{ W_mu(tau) = 1, \quad W_sigma2(tau) = z_tau / (2 * sigma). }
#'
#' @param taus Numeric vector of quantile levels (length T).
#' @param q_tau Numeric vector of baseline quantiles at taus (length T, type=8 recommended).
#' @param mu Baseline mean of Y.
#' @param sd Baseline standard deviation of Y.
#'
#' @return A T x 2 numeric matrix W with columns:
#'   \itemize{
#'     \item \code{"beta_mu"} effect of SNP on the mean
#'     \item \code{"beta_sigma2"} effect of SNP on the variance
#'   }
#' @examples
#' taus <- seq(0.1, 0.9, by = 0.2)
#' y <- rnorm(2000, mean = 2, sd = 1.5)
#' q_tau <- as.numeric(quantile(y, taus, type = 8))
#' W <- make_weight_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))
#' W
#' @export
make_weight_vqtl <- function(
  taus,
  q_tau,
  mu,
  sd
) {
  stopifnot(length(taus) == length(q_tau))
  if (sd <= 0) stop("sd must be positive.")

  z_tau <- qnorm(taus)  # baseline standard normal quantiles

  W_mu     <- rep(1, length(taus))
  W_sigma2 <- z_tau / (2 * sd)

  W <- cbind(W_mu, W_sigma2)
  colnames(W) <- c("beta_mu", "beta_sigma2")
  rownames(W) <- paste0("tau", taus)
  storage.mode(W) <- "double"
  W
}


#' @title
#' RIF quantile GWAS with generic tau-to-parameter mapping
#'
#' @details
#' Fast pipeline to: Build RIF(Y; tau) at user-specified \code{tau} via KDE of f_Y(q_tau).For each SNP, obtain OLS tau-slopes b(tau) and naive per-tau SEs. Map b(tau) to user-chosen parameters via W using OLS pseudoinverse 
#'  A = (W'W)^{-1} W' (no GLS). Compute delta-method SEs for the mapped parameters, using either a diagonal approximation or a plugin tau-correlation.
#' 
#'
#' Covariates: by default, only the SNP is residualized on \code{C}.
#' If \code{residualize_Y=TRUE}, both \code{Y} and each SNP are residualized on \code{C}.
#'
#' @param Y Numeric vector (length N), outcome.
#' @param G Numeric matrix N x P of SNP dosages (or similar; QC/I-O not performed here).
#' @param taus Numeric vector of quantile levels (default \code{seq(0.10,0.90,0.05)}).
#' @param C Optional N x K covariate matrix.
#' @param residualize_Y Logical; if TRUE, regress \code{Y} on \code{C} before quantiles (and SNPs on \code{C} too).
#' @param density_floor Floor for estimated f_Y(q_tau).
#' @param transform One of \code{"custom_W"} or \code{"two_normal"}.
#' @param transform_args List of arguments passed to the chosen transform:
#'   \itemize{
#'     \item For \code{transform="custom_W"}: provide \code{W} (T x K).
#'     \item For \code{transform="two_normal"}: provide \code{p1, mu1, sd1, mu2, sd2, include_membership} (logical).
#'   }
#' @param se_mode \code{"diagonal"} (fast, default) or \code{"plugin_cor"} for delta SEs.
#' @param plugin_R Optional T x T correlation matrix between tau-slopes; if NULL and \code{se_mode="plugin_cor"},
#'   it is estimated as \code{cor(t(Q_slope))} and projected to near-PD.
#' @param benchmark Logical; if TRUE, include timing information.
#' @param verbose Logical; toggle brief progress messages.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{taus, q_tau, fhat_tau}
#'   \item \code{Q_slope} (T x P) and \code{SE_tau} (T x P)
#'   \item \code{transform, transform_args, W} and \code{A} (pseudoinverse)
#'   \item \code{params} (K x P) and \code{SE_params} (K x P), rownames from \code{W}
#'   \item If \code{se_mode="plugin_cor"}, the \code{R_tau} used
#'   \item If \code{benchmark=TRUE}, \code{timing} sublist
#' }
#'
#' @details
#' Methodology: RIF at tau uses q_tau and f_Y(q_tau) estimated via Gaussian KDE
#' with \code{bw.nrd0}. Per-SNP tau-slopes use closed-form OLS with per-tau residual variance for naive SEs.
#' The tau-to-parameter mapping uses the OLS pseudoinverse of W. Delta-method SEs are produced either
#' under a diagonal variance (independent tau) or using a plugin tau-correlation.
#'
#' @examples
#' set.seed(42)
#' N <- 5000; P <- 50
#' taus <- seq(0.10, 0.90, by = 0.05)
#' maf <- runif(P, 0.05, 0.5)
#' G <- matrix(rbinom(N * P, 2, rep(maf, each = N)), N, P)
#' # Mixture phenotype
#' z  <- rbinom(N, 1, 0.5)
#' Y  <- rnorm(N, mean = ifelse(z==1, 3.0, 1.2), sd = ifelse(z==1, 0.7, 0.45))
#'
#' # Run with two_normal transform
#' fit <- qgwas_rif(
#'   Y, G, taus = taus,
#'   transform = "two_normal",
#'   transform_args = list(p1 = 0.5, mu1 = 1.2, sd1 = 0.45, mu2 = 3.0, sd2 = 0.7, include_membership = TRUE),
#'   se_mode = "diagonal", benchmark = FALSE, verbose = FALSE
#' )
#'
#' # The same using custom_W
#' Wsame <- make_weights_normal_mixture(taus, fit$q_tau, 0.5, 1.2, 0.45, 3.0, 0.7, include_membership = TRUE)
#' fit2  <- qgwas_rif(
#'   Y, G, taus = taus,
#'   transform = "custom_W",
#'   transform_args = list(W = Wsame),
#'   se_mode = "diagonal", benchmark = FALSE, verbose = FALSE
#' )
#' # Estimates match up to numerical tolerance:
#' max(abs(fit$params - fit2$params))
#'
#' @export
qgwas_rif <- function(
  Y, G,
  taus = seq(0.10, 0.90, 0.05),
  C = NULL,
  residualize_Y = FALSE,
  density_floor = 1e-8,
  transform = c("custom_W", "two_normal"),
  transform_args = list(),
  se_mode = c("diagonal", "plugin_cor"),
  plugin_R = NULL,
  benchmark = TRUE,
  verbose = TRUE
) {
  transform <- match.arg(transform)
  se_mode   <- match.arg(se_mode)

  Y <- as.numeric(Y)
  G <- as.matrix(G)
  N <- length(Y); P <- ncol(G); Tt <- length(taus)
  if (nrow(G) != N) stop("nrow(G) must equal length(Y).")

  t_all0 <- proc.time()[["elapsed"]]

  # 1) Residualize Y (optional) and build RIF
  if (!is.null(C) && isTRUE(residualize_Y)) {
    if (verbose) message("Residualizing Y on C...")
    Y_star <- .residualize_on_C(Y, C)
  } else {
    Y_star <- Y
  }

  if (verbose) message("Building RIF matrix...")
  t0 <- proc.time()[["elapsed"]]
  rif <- .build_rif_matrix(Y_star, taus, density_floor = density_floor)
  RIF     <- rif$R              # N x T
  q_tau   <- rif$q_tau
  fhat    <- rif$fhat_tau
  R_mean  <- colMeans(RIF)
  t1 <- proc.time()[["elapsed"]]
  time_rif <- t1 - t0

  # 2) Per-SNP tau-slopes + naive per-tau SEs
  if (verbose) message("Computing per-SNP tau-slopes and naive SEs...")
  t2 <- proc.time()[["elapsed"]]
  Q_slope <- matrix(NA_real_, nrow = Tt, ncol = P,
                    dimnames = list(paste0("tau", taus), colnames(G)))
  SE_tau  <- matrix(NA_real_, nrow = Tt, ncol = P,
                    dimnames = list(paste0("tau", taus), colnames(G)))

  # Precompute a column-replicated R_mean
  term1_base <- NULL # built per-SNP (a depends on SNP)
  for (j in seq_len(P)) {
    x <- G[, j]

    if (!is.null(C)) {
      # Always residualize SNP on C; residualize Y only if requested
      x <- .residualize_on_C(x, C)
    }

    x_c <- x - mean(x)
    Sxx <- sum(x_c * x_c)
    if (!(is.finite(Sxx) && Sxx > 0)) next

    # Slopes across tau
    b <- as.numeric(crossprod(x_c, RIF) / Sxx)

    # Intercepts a(tau) = mean(RIF_tau) - b(tau) * mean(x)
    a <- R_mean - b * mean(x)

    # Residual matrix E = RIF - (a + b * x)
    term1 <- matrix(rep(a, each = N), nrow = N, ncol = Tt)
    term2 <- outer(x, b)
    E     <- RIF - term1 - term2

    SSE    <- colSums(E * E)
    sigma2 <- SSE / max(N - 2L, 1L)
    se_b   <- sqrt(sigma2 / Sxx)

    Q_slope[, j] <- b
    SE_tau[,  j] <- se_b
  }
  t3 <- proc.time()[["elapsed"]]
  time_gwas <- t3 - t2

  # 3) Choose/construct W and its OLS pseudoinverse A
  if (transform == "custom_W") {
    if (is.null(transform_args$W)) stop("For transform='custom_W', provide transform_args$W (T x K).")
    W <- as.matrix(transform_args$W)
    if (!all(dim(W)[1] == Tt)) stop("nrow(W) must equal length(taus).")
  } else {
    # two_normal
    req <- c("p1","mu1","sd1","mu2","sd2","include_membership")
    missing_args <- setdiff(req, names(transform_args))
    if (length(missing_args))
      stop("Missing transform_args for two_normal: ", paste(missing_args, collapse = ", "))
    W <- make_weights_normal_mixture(
      taus = taus, q_tau = q_tau,
      p1 = transform_args$p1,
      mu1 = transform_args$mu1, sd1 = transform_args$sd1,
      mu2 = transform_args$mu2, sd2 = transform_args$sd2,
      include_membership = isTRUE(transform_args$include_membership)
    )
  }
  K <- ncol(W)
  A <- solve(crossprod(W), t(W))  # (K x T)

  # 4) Map tau-slopes -> parameters
  # params = A %*% Q_slope  (K x P)
  params <- A %*% Q_slope
  rownames(params) <- colnames(W)
  colnames(params) <- colnames(G)

  # 5) Delta-method SEs for mapped parameters
  if (se_mode == "diagonal") {
    # Vectorized: Var(theta_i) rowwise = A^2 %*% Var(b_tau)diag
    A2 <- A * A                       # K x T
    Vb <- SE_tau * SE_tau             # T x P
    SE_params <- sqrt(A2 %*% Vb)      # K x P
    R_used <- NULL
  } else {
    # plugin tau-correlation
    if (is.null(plugin_R)) {
      R_tau <- stats::cor(t(Q_slope), use = "pairwise.complete.obs")
    } else {
      R_tau <- as.matrix(plugin_R)
    }
    # Near-PD repair
    R_tau <- .nearPD_eig(R_tau, eps = 1e-8)

    SE_params <- matrix(NA_real_, nrow = K, ncol = P,
                        dimnames = list(colnames(W), colnames(G)))
    for (j in seq_len(P)) {
      d  <- SE_tau[, j]
      Dj <- diag(d, nrow = Tt, ncol = Tt)
      Sigma_j <- Dj %*% R_tau %*% Dj
      Vtheta  <- A %*% Sigma_j %*% t(A)
      # numerical guard
      Vtheta  <- .nearPD_eig(Vtheta, eps = 1e-12)
      SE_params[, j] <- sqrt(pmax(diag(Vtheta), 0))
    }
    R_used <- R_tau
  }

  t_all1 <- proc.time()[["elapsed"]]
  timing <- list(
    rif_build_sec = time_rif,
    gwas_loop_sec = time_gwas,
    total_sec     = t_all1 - t_all0
  )

  out <- list(
    taus       = taus,
    q_tau      = q_tau,
    fhat_tau   = fhat,
    Q_slope    = Q_slope,
    SE_tau     = SE_tau,
    transform  = transform,
    transform_args = transform_args,
    W          = W,
    A          = A,
    params     = params,
    SE_params  = SE_params
  )
  if (se_mode == "plugin_cor") out$R_tau <- R_used
  if (benchmark) out$timing <- timing
  out
}
