# ============================================================
# Proof 3: Stage 2 SE Mode Comparison & Mixture Recovery
#
# This is the third validation step in the sequence:
# 1) OLS jackknife vs analytic OLS (stage1)
# 2) RIF/quantile jackknife vs analytic RIF OLS (stage1)
# 3) Stage2 SE modes vs RIF regression (stage2)
# ============================================================

library(fungwasStage2)
library(ggplot2)
library(patchwork)
library(jsonlite)

set.seed(12345)

cat("\n======================================================\n")
cat("Running Stage 2 Proof via Stage 1 Simulation\n")
cat("======================================================\n")

# -------------------------
# 1. Run Stage 1 Simulation (Python)
# -------------------------
taus <- seq(0.1, 0.9, 0.1)
out_prefix <- tempfile("stage1_sim_")
py_script <- file.path("docs", "examples", "generate_stage1_sim.py")

cat("Simulating data + running Stage 1 (Python)...\n")
system2(
  "python",
  args = c(
    py_script,
    "--out-prefix", out_prefix,
    "--n-samples", "300000",
    "--n-snps", "100",
    "--n-covars", "2",
    "--seed", "12345",
    "--taus", paste(taus, collapse = ","),
    "--n-blocks", "25",
    "--n-threads", "4"
  ),
  stdout = TRUE
)

stage1_file <- paste0(out_prefix, ".stage1.tsv.gz")
cov_file <- paste0(out_prefix, ".cov.gz")
meta_file <- paste0(out_prefix, ".meta.json")

meta <- jsonlite::read_json(meta_file, simplifyVector = TRUE)

# -------------------------
# 2. Stage 2 with ALL SE MODES
# -------------------------
q_tau <- as.numeric(meta$q_tau)
W <- make_weights_normal_mixture(
  taus, q_tau,
  p1 = meta$p1, mu1 = meta$mu1, sd1 = meta$sd1,
  mu2 = meta$mu2, sd2 = meta$sd2,
  include_membership = TRUE
)

cat("Running Stage 2 (SE Mode: diagonal)...\n")
fit_diag <- param_gwas_from_file(stage1_file, W = W, se_mode = "diagonal")

cat("Running Stage 2 (SE Mode: plugin_cor)...\n")
fit_plugin <- param_gwas_from_file(stage1_file, W = W, se_mode = "plugin_cor")

cat("Running Stage 2 (SE Mode: dwls)...\n")
fit_dwls <- param_gwas_from_file(stage1_file, W = W, se_mode = "dwls")

cat("Running Stage 2 (SE Mode: tau_cov - GOLD STANDARD)...\n")
fit_taucov <- param_gwas_from_file(stage1_file, W = W, se_mode = "tau_cov", cov_file = cov_file)

# -------------------------
# 3. Compare SEs
# -------------------------
extract_ses <- function(fit, mode_name) {
  se_cols <- grep("_se$", names(fit), value = TRUE)
  se <- as.vector(as.matrix(fit[, se_cols, with = FALSE]))
  data.frame(SE = se, Mode = mode_name)
}

df_se <- rbind(
  extract_ses(fit_diag, "Diagonal"),
  extract_ses(fit_plugin, "Plugin"),
  extract_ses(fit_dwls, "DWLS"),
  extract_ses(fit_taucov, "TauCov (JK)")
)

se_ref <- as.vector(as.matrix(fit_taucov[, grep("_se$", names(fit_taucov), value = TRUE), with = FALSE]))
valid_se <- is.finite(se_ref) & se_ref > 0
ratio_diag <- as.vector(as.matrix(fit_diag[, grep("_se$", names(fit_diag), value = TRUE), with = FALSE]))[valid_se] / se_ref[valid_se]
ratio_plugin <- as.vector(as.matrix(fit_plugin[, grep("_se$", names(fit_plugin), value = TRUE), with = FALSE]))[valid_se] / se_ref[valid_se]
ratio_dwls <- as.vector(as.matrix(fit_dwls[, grep("_se$", names(fit_dwls), value = TRUE), with = FALSE]))[valid_se] / se_ref[valid_se]

cat("\n--- SE Comparisons (Ratio to TauCov) ---\n")
print_stats <- function(name, r) {
  cat(sprintf("%s / TauCov: Median=%.3f, Mean=%.3f\n", name, median(r, na.rm=TRUE), mean(r, na.rm=TRUE)))
}
print_stats("Diagonal", ratio_diag)
print_stats("Plugin", ratio_plugin)
print_stats("DWLS", ratio_dwls)

p1 <- ggplot(data.frame(Ref = se_ref, Est = as.vector(as.matrix(fit_diag[, grep("_se$", names(fit_diag), value = TRUE), with = FALSE]))),
             aes(x = Ref, y = Est)) +
  geom_point(alpha = 0.5) + geom_abline(color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  labs(title = "Diagonal vs TauCov (JK)", x = "TauCov SE", y = "Diagonal SE")

p2 <- ggplot(data.frame(Ref = se_ref, Est = as.vector(as.matrix(fit_plugin[, grep("_se$", names(fit_plugin), value = TRUE), with = FALSE]))),
             aes(x = Ref, y = Est)) +
  geom_point(alpha = 0.5) + geom_abline(color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  labs(title = "Plugin vs TauCov (JK)", x = "TauCov SE", y = "Plugin SE")

# -------------------------
# 4. Effect Recovery (TauCov)
# -------------------------
cat("\n--- Effect Recovery (TauCov Mode) ---\n")
param_names <- c("gamma", "beta_1", "beta_2")
true_map <- list(
  gamma = unlist(meta$beta_gamma_true),
  beta_1 = unlist(meta$beta_1_true),
  beta_2 = unlist(meta$beta_2_true)
)

for (nm in param_names) {
  est <- fit_taucov[[nm]]
  if (!is.null(est)) {
    r_val <- cor(est, true_map[[nm]])
    cat(sprintf("%s Correlation: %.4f\n", nm, r_val))
  }
}

p3 <- ggplot(data.frame(True = true_map$gamma, Est = fit_taucov$gamma), aes(x = True, y = Est)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_abline(color = "red") +
  labs(title = "Gamma Recovery", x = "True Beta", y = "Estimated Beta")

p4 <- ggplot(data.frame(True = true_map$beta_1, Est = fit_taucov$beta_1), aes(x = True, y = Est)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_abline(color = "red") +
  labs(title = "Beta_1 Recovery", x = "True Beta", y = "Estimated Beta")

final_plot <- (p1 + p2) / (p3 + p4)
ggsave("proof_stage2_results.png", final_plot, width = 10, height = 8)
cat("\nSaved comparison plot to proof_stage2_results.png\n")
