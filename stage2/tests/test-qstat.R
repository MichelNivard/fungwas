# ============================================================
# Test: Model comparison using Q statistic
# ============================================================

devtools::document()
library(fungwas)

set.seed(123)

# -------------------------
# 1. Simulate genotypes
# -------------------------
N <- 15000
P <- 500
maf <- runif(P, 0.05, 0.5)
G <- matrix(rbinom(N * P, 2, rep(maf, each = N)), nrow = N, ncol = P)
colnames(G) <- paste0("SNP", seq_len(P))

# -------------------------
# 2. Define true effects
# -------------------------
beta_mu_true     <- rnorm(P) / 40
beta_sigma2_true <- rnorm(P) / 40

# -------------------------
# 3. Simulate phenotype
# -------------------------
mu0 <- 4.0
sd0 <- 0.6

mu <- mu0 + G %*% beta_mu_true
sigma2 <- sd0^2 + G %*% beta_sigma2_true
sigma2[sigma2 <= 0] <- 1e-4  # guard against negative variances

Y <- rnorm(N, mean = mu, sd = sqrt(sigma2))

# -------------------------
# 4. Stage 1: run quantile GWAS
# -------------------------
taus <- seq(0.05, 0.95, 0.05)

stage1 <- quantile_gwas(
  Y, G,
  taus = taus,
  benchmark = FALSE,
  verbose = FALSE
)

# -------------------------
# 5. Stage 2: map to mean/variance parameters
# -------------------------
W_var <- make_weights_vqtl(taus, stage1$q_tau, mu = mean(Y), sd = sd(Y))

fitcorrect <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_var),
  se_mode = "dwls"
)

# -------------------------
# 4. Stage 2: wrong model (log-normal)
# -------------------------

# Estimate weights
fitlog <- fitdist(Y, "lnorm")


q_tau <- as.numeric(quantile(Y, taus, type = 8))

dist_cdf <- function(y, params) plnorm(y, meanlog = params$meanlog, sdlog = params$sdlog)
dist_pdf <- function(y, params) dlnorm(y, meanlog = params$meanlog, sdlog = params$sdlog)

###  Notice I sign fliped the derivatives...

grad_meanlog <- function(y, params) {
  z <- (log(y) - params$meanlog) / params$sdlog
  dnorm(z) / params$sdlog
}

W_analytic <- make_weights_generic(
  taus, q_tau, dist_cdf, dist_pdf,
  params = list(meanlog = fitlog$estimate[1], sdlog = fitlog$estimate[2]),
  grad_funcs = list(beta_meanlog = grad_meanlog)
)

grad_fd_meanlog <- make_fd_grad("meanlog", dist_cdf)
grad_fd_sdlog   <- make_fd_grad("sdlog", dist_cdf)


W_fd <- make_weights_generic(
  taus, q_tau, dist_cdf, dist_pdf,
  params = list(meanlog = fitlog$estimate[1], sdlog = fitlog$estimate[2]),
  grad_funcs = list(beta_meanlog = grad_fd_meanlog)
)


fit_wrong <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_analytic),
  se_mode =  "dwls"
)

# -------------------------
# 5. Inspect Q statistics
# -------------------------
cat("\n--- Model Comparison ---\n")
cat("Correct model Q (per SNP):\n")
print(summary(fit_correct$Q))

cat("\nWrong model Q (per SNP):\n")
print(summary(fit_wrong$Q))

cat("\nAverage Q difference (wrong - correct):\n")

plot(fit_correct$Q-fit_wrong$Q)
