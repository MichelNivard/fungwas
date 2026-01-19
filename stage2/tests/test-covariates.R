# ============================================================
# Test script: Confounding via shared latent factor
# First 60 SNPs have allele frequencies + phenotype influenced by same U
# ============================================================

devtools::document()
library(fungwasStage2)
library(ggplot2)

set.seed(42)

# -------------------------
# 1. Simulation parameters
# -------------------------
N <- 100000
P <- 500
maf <- runif(P, 0.05, 0.5)  # baseline allele frequencies

# latent confounder (e.g. ancestry or environment)
U <- rnorm(N)

# -------------------------
# 2. Simulate genotypes
# -------------------------
G <- matrix(NA, nrow = N, ncol = P)

# First 60 SNPs: allele frequency depends on U (confounded)
for (j in 1:60) {
  p_j <- plogis(qlogis(maf[j]) + 0.3 * scale(U))  # logit shift of AF by U
  G[, j] <- rbinom(N, 2, p_j)
}

# Remaining SNPs: independent of U
for (j in 61:P) {
  G[, j] <- rbinom(N, 2, maf[j])
}

colnames(G) <- paste0("SNP", seq_len(P))

# -------------------------
# 3. True SNP effects on Y
# -------------------------
beta_mu_true     <- rnorm(P) / 40
beta_sigma2_true <- rnorm(P) / 40

# -------------------------
# 4. Phenotype: affected by both SNPs and confounder U
# -------------------------
mu0 <- 7.0
sd0 <- 1.25

mu <- mu0 + G %*% beta_mu_true + 0.5 * U        # U → mean
sigma2 <- sd0^2 + G %*% beta_sigma2_true + 0.2 * U  # U → variance
sigma2[sigma2 <= 0] <- 1e-4

Y <- rnorm(N, mean = mu, sd = sqrt(sigma2))

# -------------------------
# 5. Stage 1: run quantile GWAS
# -------------------------
taus <- seq(0.1, 0.9, 0.1)

# (a) without covariates (naive GWAS)
stage1_naive <- quantile_gwas(
  Y, G,
  taus = taus,
  benchmark = FALSE,
  verbose = FALSE
)

# (b) with U as covariate (adjusted GWAS)
stage1_adj <- quantile_gwas(
  Y, G,
  C = cbind(U),
  taus = taus,
  benchmark = FALSE,
  verbose = FALSE,
  covar_mode = "residualize"
)

# -------------------------
# 6. Map to mean/variance parameters
# -------------------------
W_var <- make_weights_vqtl(taus, stage1_naive$q_tau, mu = mean(Y), sd = sd(Y))

stage2_naive <- param_gwas(
  stage1_naive,
  transform = "custom_W",
  transform_args = list(W = W_var),
  se_mode = "plugin_cor"
)

stage2_adj <- param_gwas(
  stage1_adj,
  transform = "custom_W",
  transform_args = list(W = W_var),
  se_mode = "plugin_cor"
)

# -------------------------
# 7. Inspect SNP effects
# -------------------------
est_naive <- t(stage2_naive$params)
est_adj   <- t(stage2_adj$params)

results <- data.frame(
  SNP              = colnames(G),
  beta_mu_true     = beta_mu_true,
  beta_mu_hat_naive = est_naive[, "beta_mu"],
  beta_mu_hat_adj   = est_adj[, "beta_mu"]
)

# -------------------------
# 8. Evaluate bias
# -------------------------
cor_naive <- cor(results$beta_mu_true, results$beta_mu_hat_naive)
cor_adj   <- cor(results$beta_mu_true, results$beta_mu_hat_adj)

cat(sprintf("Correlation true vs est (mean effects):\n"))
cat(sprintf("  Naive (no covariates): %.3f\n", cor_naive))
cat(sprintf("  Adjusted (with U):     %.3f\n", cor_adj))

# -------------------------
# 9. Plots: highlight first 60 SNPs (confounded)
# -------------------------
results$group <- ifelse(seq_len(nrow(results)) <= 60, "Confounded", "Other")

p_mu <- ggplot(results, aes(x = beta_mu_true, y = beta_mu_hat_naive, color = group)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Confounded" = "blue", "Other" = "steelblue")) +
  labs(title = "Naive GWAS: Mean effects (confounded SNPs in blue)",
       x = "True beta_mu", y = "Estimated beta_mu") +
  theme_minimal()

p_mu_adj <- ggplot(results, aes(x = beta_mu_true, y = beta_mu_hat_adj, color = group)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Confounded" = "blue", "Other" = "steelblue")) +
  labs(title = "Adjusted GWAS: Mean effects (confounded SNPs in blue)",
       x = "True beta_mu", y = "Estimated beta_mu") +
  theme_minimal()

print(p_mu)
print(p_mu_adj)


# -------------------------
# Variance effects (naive vs adjusted)
# -------------------------

results$beta_sigma2_true      <- beta_sigma2_true
results$beta_sigma2_hat_naive <- t(stage2_naive$params)[, "beta_sigma2"]
results$beta_sigma2_hat_adj   <- t(stage2_adj$params)[, "beta_sigma2"]

# Naive variance effect plot
p_sigma2_naive <- ggplot(results, aes(x = beta_sigma2_true, 
                                      y = beta_sigma2_hat_naive, 
                                      color = group)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Confounded" = "blue", "Other" = "darkorange")) +
  labs(title = "Naive GWAS: Variance effects (confounded SNPs in blue)",
       x = "True beta_sigma2", y = "Estimated beta_sigma2") +
  theme_minimal()

# Adjusted variance effect plot
p_sigma2_adj <- ggplot(results, aes(x = beta_sigma2_true, 
                                    y = beta_sigma2_hat_adj, 
                                    color = group)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Confounded" = "blue", "Other" = "darkorange")) +
  labs(title = "Adjusted GWAS: Variance effects (confounded SNPs in blue)",
       x = "True beta_sigma2", y = "Estimated beta_sigma2") +
  theme_minimal()

print(p_sigma2_naive)
print(p_sigma2_adj)
