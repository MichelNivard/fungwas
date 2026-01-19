# ============================================================
# Test script for two-component Normal mixture with mean, variance, and membership effects
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
beta_gamma_true  <- rnorm(P) / 40    # effects on membership (logit p1)
beta_mu1_true    <- rnorm(P) / 40    # effects on mean of component 1
beta_mu2_true    <- rnorm(P) / 40    # effects on mean of component 2
beta_sigma1_true <- rnorm(P) / 60    # effects on sd of component 1
beta_sigma2_true <- rnorm(P) / 60    # effects on sd of component 2

# -------------------------
# 3. Baseline mixture parameters
# -------------------------
p1  <- 0.4
mu1 <- 1.0; sd1 <- 0.5
mu2 <- 3.0; sd2 <- 1.0

# -------------------------
# 4. Simulate phenotype
# -------------------------
gamma <- rep(qlogis(p1), N) + G %*% beta_gamma_true
p1_i  <- plogis(gamma)  # individual-specific class prob

# draw latent class
class <- rbinom(N, 1, p1_i)

# class-specific means and sds
mu_i <- ifelse(class == 1,
               mu1 + G %*% beta_mu1_true,
               mu2 + G %*% beta_mu2_true)

sd_i <- ifelse(class == 1,
               pmax(sd1 + G %*% beta_sigma1_true, 1e-4),
               pmax(sd2 + G %*% beta_sigma2_true, 1e-4))

Y <- rnorm(N, mean = mu_i, sd = sd_i)

# -------------------------
# 5. Stage 1: run quantile GWAS
# -------------------------
taus <- seq(0.04, 0.96, 0.02)

stage1 <- quantile_gwas(
  Y, G,
  taus = taus,
  benchmark = FALSE,
  verbose = FALSE
)

# -------------------------
# 6. Stage 2: map tau slopes into mixture-vQTL parameters
# -------------------------
W_mix_vqtl <- make_weights_normal_mixture_vqtl(
  taus, stage1$q_tau,
  p1 = p1, mu1 = mu1, sd1 = sd1,
  mu2 = mu2, sd2 = sd2
)

stage2 <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_mix_vqtl),
  se_mode = "plugin_cor"
)

# -------------------------
# 7. Collect results
# -------------------------
est <- t(stage2$params)
se  <- t(stage2$SE_params)

results <- data.frame(
  SNP               = colnames(G),
  beta_gamma_true   = beta_gamma_true,
  beta_mu1_true     = beta_mu1_true,
  beta_mu2_true     = beta_mu2_true,
  beta_sigma1_true  = beta_sigma1_true,
  beta_sigma2_true  = beta_sigma2_true,
  beta_gamma_hat    = est[, "gamma"],
  beta_mu1_hat      = est[, "beta_mu1"],
  beta_mu2_hat      = est[, "beta_mu2"],
  beta_sigma1_hat   = est[, "beta_sigma1"],
  beta_sigma2_hat   = est[, "beta_sigma2"],
  se_gamma          = se[, "gamma"],
  se_mu1            = se[, "beta_mu1"],
  se_mu2            = se[, "beta_mu2"],
  se_sigma1         = se[, "beta_sigma1"],
  se_sigma2         = se[, "beta_sigma2"]
)

# -------------------------
# 8. Evaluate recovery
# -------------------------
cor_gamma  <- cor(results$beta_gamma_true,  results$beta_gamma_hat)
cor_mu1    <- cor(results$beta_mu1_true,    results$beta_mu1_hat)
cor_mu2    <- cor(results$beta_mu2_true,    results$beta_mu2_hat)
cor_sigma1 <- cor(results$beta_sigma1_true, results$beta_sigma1_hat)
cor_sigma2 <- cor(results$beta_sigma2_true, results$beta_sigma2_hat)

cat(sprintf("Correlation (true vs est) gamma:   %.3f\n", cor_gamma))
cat(sprintf("Correlation (true vs est) mu1:     %.3f\n", cor_mu1))
cat(sprintf("Correlation (true vs est) mu2:     %.3f\n", cor_mu2))
cat(sprintf("Correlation (true vs est) sigma1:  %.3f\n", cor_sigma1))
cat(sprintf("Correlation (true vs est) sigma2:  %.3f\n", cor_sigma2))

# -------------------------
# 9. Plots
# -------------------------
plot_effects <- function(df, true, hat, title, col) {
  ggplot(df, aes_string(x = true, y = hat)) +
    geom_point(color = col, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = title,
         x = paste("True", true),
         y = paste("Estimated", hat)) +
    theme_minimal()
}

print(plot_effects(results, "beta_gamma_true",  "beta_gamma_hat",  "Membership effects", "purple"))
print(plot_effects(results, "beta_mu1_true",    "beta_mu1_hat",    "Mean effects (comp 1)", "steelblue"))
print(plot_effects(results, "beta_mu2_true",    "beta_mu2_hat",    "Mean effects (comp 2)", "darkorange"))
print(plot_effects(results, "beta_sigma1_true", "beta_sigma1_hat", "SD effects (comp 1)", "forestgreen"))
print(plot_effects(results, "beta_sigma2_true", "beta_sigma2_hat", "SD effects (comp 2)", "brown"))
