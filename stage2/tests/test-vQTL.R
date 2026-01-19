# ============================================================
# Test script for vQTL detection with fgwas
# ============================================================

devtools::document()
library(fungwas)
library(ggplot2)

set.seed(123)

# -------------------------
# 1. Simulate genotypes
# -------------------------
N <- 150000
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
mu0 <- 2.0
sd0 <- 1.0

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

stage2 <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_var),
  se_mode = "plugin_cor"
)

# -------------------------
# 6. Inspect results
# -------------------------
est <- t(stage2$params)    # P x 2, columns = beta_mu, beta_sigma2
se  <- t(stage2$SE_params)

results <- data.frame(
  SNP              = colnames(G),
  beta_mu_true     = beta_mu_true,
  beta_sigma2_true = beta_sigma2_true,
  beta_mu_hat      = est[, "beta_mu"],
  beta_sigma2_hat  = est[, "beta_sigma2"],
  se_mu            = se[, "beta_mu"],
  se_sigma2        = se[, "beta_sigma2"]
)

# -------------------------
# 7. Evaluate recovery
# -------------------------
cor_mu     <- cor(results$beta_mu_true, results$beta_mu_hat)
cor_sigma2 <- cor(results$beta_sigma2_true, results$beta_sigma2_hat)

cat(sprintf("Correlation (true vs est) mean effects:     %.3f\n", cor_mu))
cat(sprintf("Correlation (true vs est) variance effects: %.3f\n", cor_sigma2))

# -------------------------
# 8. Plots
# -------------------------
p_mu <- ggplot(results, aes(x = beta_mu_true, y = beta_mu_hat)) +
  geom_point(color = "steelblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Mean effects", x = "True beta_mu", y = "Estimated beta_mu") +
  theme_minimal()

p_sigma2 <- ggplot(results, aes(x = beta_sigma2_true, y = beta_sigma2_hat)) +
  geom_point(color = "darkorange", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Variance effects", x = "True beta_sigma2", y = "Estimated beta_sigma2") +
  theme_minimal()

print(p_mu)
print(p_sigma2)
