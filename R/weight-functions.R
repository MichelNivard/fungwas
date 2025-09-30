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

#' Build tau-by-K weights for a two-component Normal mixture with mean, variance, and membership effects
#'
#' Constructs a weight matrix \code{W} that maps RIF tau-slopes \eqn{b(tau)}
#' into SNP effects on the parameters of a two-component Normal mixture:
#' \deqn{Y ~ p1 * N(mu1, sd1^2) + (1 - p1) * N(mu2, sd2^2).}
#'
#' Supported SNP effect targets:
#' \itemize{
#'   \item \code{gamma}: log-odds of membership in component 1 (class proportion).
#'   \item \code{beta_mu1}, \code{beta_mu2}: effects on the component means.
#'   \item \code{beta_sigma1}, \code{beta_sigma2}: effects on the component standard deviations.
#' }
#'
#' @details
#' For baseline quantiles \eqn{q_tau} and mixture pdf
#' \deqn{f(q) = p1 f1(q) + (1 - p1) f2(q),}
#' where \eqn{f1, f2} are the component pdfs and \eqn{F1, F2} their CDFs, the weights are:
#'
#' - Membership:
#'   \deqn{W_gamma(tau) = - (p1 * (1 - p1) / f(q_tau)) * (F1(q_tau) - F2(q_tau))}
#'
#' - Means:
#'   \deqn{W_mu1(tau) = (p1 * f1(q_tau)) / f(q_tau), \quad
#'         W_mu2(tau) = ((1 - p1) * f2(q_tau)) / f(q_tau)}
#'
#' - Standard deviations:
#'   \deqn{W_sigma1(tau) = (p1 * f1(q_tau) * (q_tau - mu1) / sd1) / f(q_tau),}
#'   \deqn{W_sigma2(tau) = ((1 - p1) * f2(q_tau) * (q_tau - mu2) / sd2) / f(q_tau).}
#'
#' These arise from the derivative of the mixture CDF with respect to each parameter,
#' mapped into the quantile effect scale via the implicit function theorem.
#'
#' @param taus Numeric vector of quantile levels (length T).
#' @param q_tau Numeric vector of baseline quantiles of Y at \code{taus} (length T, \code{type=8} recommended).
#' @param p1 Proportion of component 1 (in (0,1)); component 2 proportion = 1 - p1.
#' @param mu1,sd1 Mean and sd of component 1.
#' @param mu2,sd2 Mean and sd of component 2.
#' @param tiny Small positive floor to stabilize divisions (mixture pdf clamped to at least \code{tiny}).
#'
#' @return A T x 5 numeric matrix \code{W} with columns:
#'   \code{c("gamma","beta_mu1","beta_mu2","beta_sigma1","beta_sigma2")}.
#'
#' @examples
#' taus <- seq(0.1, 0.9, 0.1)
#' set.seed(1)
#' y <- c(rnorm(500, 1.2, 0.6), rnorm(500, 3.0, 0.9))
#' q_tau <- as.numeric(quantile(y, taus, type = 8))
#' W <- make_weights_normal_mixture_vqtl(
#'   taus, q_tau,
#'   p1 = 0.5, mu1 = 1.2, sd1 = 0.6,
#'   mu2 = 3.0, sd2 = 0.9
#' )
#' head(W)
#' @export
make_weights_normal_mixture_vqtl <- function(
  taus,
  q_tau,
  p1,
  mu1, sd1,
  mu2, sd2,
  tiny = 1e-12
) {
  stopifnot(length(taus) == length(q_tau))

  # Component pdfs and cdfs
  f1 <- dnorm(q_tau, mean = mu1, sd = sd1)
  F1 <- pnorm(q_tau, mean = mu1, sd = sd1)
  f2 <- dnorm(q_tau, mean = mu2, sd = sd2)
  F2 <- pnorm(q_tau, mean = mu2, sd = sd2)

  # Mixture density (guarded)
  f <- pmax(p1 * f1 + (1 - p1) * f2, tiny)

  # Weights
  W_gamma    <- -(p1 * (1 - p1) / f) * (F1 - F2)
  W_mu1      <- (p1 * f1) / f
  W_mu2      <- ((1 - p1) * f2) / f
  W_sigma1   <- (p1 * f1 * (q_tau - mu1) / sd1) / f
  W_sigma2   <- ((1 - p1) * f2 * (q_tau - mu2) / sd2) / f

  W <- cbind(W_gamma, W_mu1, W_mu2, W_sigma1, W_sigma2)
  colnames(W) <- c("gamma","beta_mu1","beta_mu2","beta_sigma1","beta_sigma2")
  rownames(W) <- paste0("tau", taus)
  storage.mode(W) <- "double"
  W
}
