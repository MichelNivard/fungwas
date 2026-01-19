# ============================================================
# Test script: Log-normal GWAS (meanlog, sdlog)
# ============================================================

devtools::document()
library(fungwas)
library(ggplot2)

set.seed(123)

# -------------------------
# 1. Simulate genotypes
# -------------------------
N <- 10000
P <- 100
maf <- runif(P, 0.05, 0.5)
G <- matrix(rbinom(N * P, 2, rep(maf, each = N)), nrow = N, ncol = P)
colnames(G) <- paste0("SNP", seq_len(P))

# -------------------------
# 2. Define true effects
# -------------------------
beta_meanlog_true <- rnorm(P) / 40
beta_sdlog_true   <- rnorm(P) / 40

# -------------------------
# 3. Simulate phenotype
# -------------------------
meanlog0 <- 0.5
sdlog0   <- 0.7

meanlog <- meanlog0 + G %*% beta_meanlog_true
sdlog   <- pmax(sdlog0 + G %*% beta_sdlog_true, 0.01)

Y <- rlnorm(N, meanlog = meanlog, sdlog = sdlog)

# -------------------------
# 4. Build weights (analytic and FD)
# -------------------------
taus <- seq(0.1, 0.9, 0.05)
q_tau <- as.numeric(quantile(Y, taus, type = 8))

dist_cdf <- function(y, params) plnorm(y, meanlog = params$meanlog, sdlog = params$sdlog)
dist_pdf <- function(y, params) dlnorm(y, meanlog = params$meanlog, sdlog = params$sdlog)

###  NOtice I sign fliped the derivatives...

grad_meanlog <- function(y, params) {
  z <- (log(y) - params$meanlog) / params$sdlog
  dnorm(z) / params$sdlog
}
grad_sdlog <- function(y, params) {
  z <- (log(y) - params$meanlog) / params$sdlog
  (z * dnorm(z)) / params$sdlog
}

W_analytic <- make_weights_generic(
  taus, q_tau, dist_cdf, dist_pdf,
  params = list(meanlog = meanlog0, sdlog = sdlog0),
  grad_funcs = list(beta_meanlog = grad_meanlog, beta_sdlog = grad_sdlog)
)

grad_fd_meanlog <- make_fd_grad("meanlog", dist_cdf)
grad_fd_sdlog   <- make_fd_grad("sdlog", dist_cdf)

W_fd <- make_weights_generic(
  taus, q_tau, dist_cdf, dist_pdf,
  params = list(meanlog = meanlog0, sdlog = sdlog0),
  grad_funcs = list(beta_meanlog = grad_fd_meanlog, beta_sdlog = grad_fd_sdlog)
)

# -------------------------
# 5. Stage 1: quantile GWAS
# -------------------------
stage1 <- quantile_gwas(
  Y, G,
  taus = taus,
  benchmark = FALSE,
  verbose = FALSE
)

# -------------------------
# 6. Stage 2: param GWAS
# -------------------------
fit_analytic <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_analytic),
  se_mode = "diagonal"
)

fit_fd <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_fd),
  se_mode = "diagonal"
)

# -------------------------
# 7. Collect results
# -------------------------
est_analytic <- t(fit_analytic$params)
est_fd       <- t(fit_fd$params)

results <- data.frame(
  SNP               = colnames(G),
  beta_meanlog_true = beta_meanlog_true,
  beta_sdlog_true   = beta_sdlog_true,
  beta_meanlog_analytic = est_analytic[, "beta_meanlog"],
  beta_sdlog_analytic   = est_analytic[, "beta_sdlog"],
  beta_meanlog_fd       = est_fd[, "beta_meanlog"],
  beta_sdlog_fd         = est_fd[, "beta_sdlog"]
)

# -------------------------
# 8. Compare recovery
# -------------------------
cat("Correlation with true effects:\n")
cat(sprintf("Meanlog (analytic): %.3f\n", cor(results$beta_meanlog_true, results$beta_meanlog_analytic)))
cat(sprintf("Sdlog   (analytic): %.3f\n", cor(results$beta_sdlog_true,   results$beta_sdlog_analytic)))
cat(sprintf("Meanlog (FD):       %.3f\n", cor(results$beta_meanlog_true, results$beta_meanlog_fd)))
cat(sprintf("Sdlog   (FD):       %.3f\n", cor(results$beta_sdlog_true,   results$beta_sdlog_fd)))

cat("\nCorrelation analytic vs FD estimates:\n")
cat(sprintf("Meanlog: %.3f\n", cor(results$beta_meanlog_analytic, results$beta_meanlog_fd)))
cat(sprintf("Sdlog:   %.3f\n", cor(results$beta_sdlog_analytic, results$beta_sdlog_fd)))

# -------------------------
# 9. Plots
# -------------------------
p1 <- ggplot(results, aes(x = beta_meanlog_true, y = beta_meanlog_analytic)) +
  geom_point(color = "steelblue", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "LogNormal meanlog: true vs estimated (analytic)") +
  theme_minimal()

p2 <- ggplot(results, aes(x = beta_sdlog_true, y = beta_sdlog_analytic)) +
  geom_point(color = "darkorange", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "LogNormal sdlog: true vs estimated (analytic)") +
  theme_minimal()

print(p1)
print(p2)
