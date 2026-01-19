# ============================================================
# Test script: Shifted + scaled LogNormal
#
# Model (per SNP j):
#   W | G_j = a(G_j) + exp(beta0 + beta1_j * G_j) * U0
#   a(G_j) = a0 + a1_j * G_j
#   U0 ~ LogNormal(muU, sdU)
#
# Equivalent distribution:
#   W | G_j ~ Shifted-LogNormal( shift = a0 + a1_j G_j,
#                                meanlog = (muU + beta0) + beta1_j G_j,
#                                sdlog = sdU )
#
# Pipeline:
#   1) quantile_gwas()
#   2) make_weights_generic() for shifted LogNormal (shift + meanlog params)
#   3) param_gwas()
# ============================================================

devtools::document()
library(fungwasStage2)
library(ggplot2)

set.seed(123)

# -------------------------
# 1) Simulate genotypes
# -------------------------
N <- 60000
P <- 400
maf <- runif(P, 0.05, 0.5)
G <- matrix(rbinom(N * P, 2, rep(maf, each = N)), nrow = N, ncol = P)
colnames(G) <- paste0("SNP", 1:P)

# -------------------------
# 2) Baseline distribution + true SNP effects
# -------------------------
a0    <- 30
muU   <- 3.0
sdU   <- 0.35
beta0 <- 0.0

# Per-SNP effects on:
#   shift a(G)  : a1_j
#   log-scale   : beta1_j   (via exp(beta0 + beta1_j * G))
a1_true    <- c(rnorm(P/2, 0, sqrt(0.02)),rep(0,P/2))   # additive shift effects
beta1_true <- c(rep(0,P/2),rnorm(P/2, 0, 1/100))   # multiplicative scale effects (on log-scale)

# Generate phenotype in a "polygenic" way
aG   <- as.numeric(a0 + G %*% a1_true)
logS <- as.numeric(beta0 + G %*% beta1_true)
S    <- exp(logS)

U0 <- rlnorm(N, meanlog = muU, sdlog = sdU)
Y  <- aG + S * U0

# Ensure no numerical issues (should already be > aG)
stopifnot(all(is.finite(Y)), all(Y > 0))

# -------------------------
# 3) Stage 1: quantile GWAS
# -------------------------
taus <- seq(0.05, 0.95, 0.02)
stage1 <- quantile_gwas(Y, G, taus = taus, benchmark = FALSE, verbose = FALSE)

# -------------------------
# 4) Stage 2: Build W for Shifted LogNormal (shift + meanlog)
# -------------------------
# CDF/PDF of shifted LogNormal:
#   F(y) = 0 for y <= a
#        = plnorm(y-a, meanlog=m, sdlog=s) for y>a
#   f(y) = dlnorm(y-a, ...) for y>a
dist_cdf <- function(y, params) {
  a <- params$shift
  out <- numeric(length(y))
  ok <- (y > a)
  out[ok] <- plnorm(y[ok] - a, meanlog = params$meanlog, sdlog = params$sdlog)
  out
}

dist_pdf <- function(y, params) {
  a <- params$shift
  out <- numeric(length(y))
  ok <- (y > a)
  out[ok] <- dlnorm(y[ok] - a, meanlog = params$meanlog, sdlog = params$sdlog)
  out
}

# Gradients of the *CDF* w.r.t. parameters (for make_weights_generic)
# 1) w.r.t shift a:
#    F(y) = F0(y-a)  => dF/da = -f0(y-a)
grad_shift <- function(y, params) {
  a <- params$shift
  out <- numeric(length(y))
  ok <- (y > a)
  out[ok] <- -dlnorm(y[ok] - a, meanlog = params$meanlog, sdlog = params$sdlog)
  out
}

# 2) w.r.t meanlog m (same as ordinary LogNormal, evaluated at y-a)
grad_meanlog <- function(y, params) {
  a <- params$shift
  out <- numeric(length(y))
  ok <- (y > a)
  yy <- y[ok] - a
  z  <- (log(yy) - params$meanlog) / params$sdlog
  out[ok] <- -dnorm(z) / params$sdlog
  out
}

# Plug-in baseline params at G=0
params0 <- list(
  shift   = a0,
  meanlog = muU + beta0,
  sdlog   = sdU
)

W <- make_weights_generic(
  taus = taus,
  q_tau = stage1$q_tau,
  dist_cdf = dist_cdf,
  dist_pdf = dist_pdf,
  params = params0,
  grad_funcs = list(beta_shift = grad_shift, beta_meanlog = grad_meanlog)
)

fit <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W),
  se_mode = "plugin_cor"
)

# -------------------------
# 5) Evaluate recovery
# -------------------------
est <- t(fit$params)  # P x 2

res <- data.frame(
  SNP = colnames(G),
  beta_shift_true   = a1_true,
  beta_meanlog_true = beta1_true,
  beta_shift_hat    = est[, "beta_shift"],
  beta_meanlog_hat  = est[, "beta_meanlog"]
)

message("cor(beta_shift):   ", round(cor(res$beta_shift_true,   res$beta_shift_hat), 3))
message("cor(beta_meanlog): ", round(cor(res$beta_meanlog_true, res$beta_meanlog_hat), 3))

p1 <- ggplot(res, aes(beta_shift_true, beta_shift_hat)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Shift (additive) effects: true vs estimated",
       x = "True a1", y = "Estimated beta_shift") +
  theme_minimal()

p2 <- ggplot(res, aes(beta_meanlog_true, beta_meanlog_hat)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Scale (multiplicative) effects via meanlog: true vs estimated",
       x = "True beta1", y = "Estimated beta_meanlog") +
  theme_minimal()

print(p1)
print(p2)
