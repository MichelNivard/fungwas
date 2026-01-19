# ============================================================
# Test script for two-component Normal mixture GWAS
# ============================================================

devtools::document()
library(fungwasStage2)
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
beta_mu1_true   <- rnorm(P) / 40
beta_mu2_true   <- rnorm(P) / 40
beta_gamma_true <- rnorm(P) / 40

# -------------------------
# 3. Baseline mixture
# -------------------------
p1  <- 0.7
mu1 <- 1.2; sd1 <- 0.6
mu2 <- 2.5; sd2 <- 0.6

# -------------------------
# 4. Simulate phenotype
# -------------------------
gamma <- rep(qlogis(p1), N) + G %*% beta_gamma_true
p1_i  <- plogis(gamma)  # individual-specific class prob

class <- rbinom(N, 1, p1_i)

mu_i <- ifelse(class == 1,
               mu1 + G %*% beta_mu1_true,
               mu2 + G %*% beta_mu2_true)

sd_i <- ifelse(class == 1, sd1, sd2)

Y <- rnorm(N, mean = mu_i, sd = sd_i)

# -------------------------
# 5. Stage 1: run quantile GWAS
# -------------------------
taus <- seq(0.1, 0.9, 0.05)

stage1 <- quantile_gwas(
  Y, G,
  taus = taus,
  benchmark = FALSE,
  verbose = FALSE
)

# -------------------------
# 6. Stage 2: map to mixture parameters
# -------------------------
stage2 <- param_gwas(
  stage1,
  transform = "two_normal",
  transform_args = list(
    p1 = p1, mu1 = mu1, sd1 = sd1,
    mu2 = mu2, sd2 = sd2,
    include_membership = TRUE
  ),
  se_mode = "plugin_cor"
)

# -------------------------
# 7. Collect results
# -------------------------
est <- t(stage2$params)   # P x 3: gamma, beta_1, beta_2
se  <- t(stage2$SE_params)

results <- data.frame(
  SNP            = colnames(G),
  beta_gamma_true= beta_gamma_true,
  beta_mu1_true  = beta_mu1_true,
  beta_mu2_true  = beta_mu2_true,
  beta_gamma_hat = est[, "gamma"],
  beta_mu1_hat   = est[, "beta_1"],
  beta_mu2_hat   = est[, "beta_2"],
  se_gamma       = se[, "gamma"],
  se_mu1         = se[, "beta_1"],
  se_mu2         = se[, "beta_2"]
)

# -------------------------
# 8. Evaluate recovery
# -------------------------
cor_gamma <- cor(results$beta_gamma_true, results$beta_gamma_hat)
cor_mu1   <- cor(results$beta_mu1_true, results$beta_mu1_hat)
cor_mu2   <- cor(results$beta_mu2_true, results$beta_mu2_hat)

cat(sprintf("Correlation (true vs est) class membership: %.3f\n", cor_gamma))
cat(sprintf("Correlation (true vs est) mean comp 1:      %.3f\n", cor_mu1))
cat(sprintf("Correlation (true vs est) mean comp 2:      %.3f\n", cor_mu2))

# -------------------------
# 9. Plots
# -------------------------

# Scale to the logistic regression residual variance
p_gamma <- ggplot(results, aes(x = beta_gamma_true, y = beta_gamma_hat/sqrt(pi^2/3))) +
  geom_point(color = "purple", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Class membership effects", x = "True gamma", y = "Estimated gamma") +
  theme_minimal()

p_mu1 <- ggplot(results, aes(x = beta_mu1_true, y = beta_mu1_hat)) +
  geom_point(color = "steelblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Mean effects (component 1)", x = "True beta_mu1", y = "Estimated beta_mu1") +
  theme_minimal()

p_mu2 <- ggplot(results, aes(x = beta_mu2_true, y = beta_mu2_hat)) +
  geom_point(color = "darkorange", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Mean effects (component 2)", x = "True beta_mu2", y = "Estimated beta_mu2") +
  theme_minimal()

print(p_gamma)
print(p_mu1)
print(p_mu2)
