# ============================================================
# Test: Model tournament AIC/BIC preference for true model
# ============================================================

library(fungwasStage2)

set.seed(2024)

P <- 500
Tt <- 7
taus <- seq(0.1, 0.9, length.out = Tt)

# Additive-only model
W_add <- matrix(1, nrow = Tt, ncol = 1)

# Additive + multiplicative model
W_add_mult <- cbind(
  intercept = 1,
  mult = seq(-1, 1, length.out = Tt)
)

theta_true <- rbind(
  intercept = rnorm(P, 0.1, 0.03),
  mult = rnorm(P, 0.2, 0.05)
)

beta_true <- W_add_mult %*% theta_true
noise <- matrix(rnorm(Tt * P, sd = 0.005), nrow = Tt)
Q_slope <- beta_true + noise

SE_tau <- matrix(0.02, nrow = Tt, ncol = P)

stage1 <- list(Q_slope = Q_slope, SE_tau = SE_tau, taus = taus)

timing <- system.time({
  res <- param_gwas_multi(
    stage1,
    W_list = list(add = W_add, add_mult = W_add_mult),
    N = 20000,
    n_threads = 2
  )
})
cat("Model tournament time (sec):", timing[["elapsed"]], "\n")

aic_add <- res[, "AIC_add"]
bic_add <- res[, "BIC_add"]
aic_add_mult <- res[, "AIC_add_mult"]
bic_add_mult <- res[, "BIC_add_mult"]

cat("Mean AIC (add):", mean(aic_add), "\n")
cat("Mean AIC (add_mult):", mean(aic_add_mult), "\n")
cat("Mean BIC (add):", mean(bic_add), "\n")
cat("Mean BIC (add_mult):", mean(bic_add_mult), "\n")

stopifnot(mean(aic_add_mult < aic_add) > 0.9)
stopifnot(mean(bic_add_mult < bic_add) > 0.9)
