# FungWas Stage 2: Parametric Mapping (R Package)

This R package implements "Stage 2" of the FungWas framework: mapping quantile GWAS results (from Stage 1) to parameters of a user-chosen structural model.

**Core idea:** FungWAS is a general framework for Quantile GWAS. You fit **your** null phenotype model externally (normal, mixture, skewed, etc.), then pass the resulting parameters (e.g., posterior weights, residuals, or model-derived weights `W`) into Stage 2. Mixture models are a common use case (e.g., Alfred's serology data), but they are **not** required.

## Key Features
-   **Flexible Mapping:** Map quantile effects to ANY parametric model where you can derive $dQ(\tau) / d\theta$.
-   **Accurate SEs:** Uses rigorous Jackknife-based covariance from Stage 1 (`tau_cov` mode) or efficient approximations.
-   **Pre-built Models:** Includes helpers for:
    -   vQTLs (`make_weights_vqtl`)
    -   Normal Mixtures (`make_weights_normal_mixture`)
    -   Log-normal shifts (`make_weights_lognormal`)

## Installation

```r
devtools::install(".")
library(fungwasStage2)
```

## Workflow: From Raw Data to GWAS

### Step 0: The Research Question
Examples of **why** you would use Quantile GWAS:

- **Example A (Standard Model):** "Does a SNP affect the variance of BMI?"
- **Example B (Mixture Model):** "Does a SNP affect the risk of being seropositive in a continuous antibody distribution?"

### Step 1: Fit Your Model (External to FungWAS)
You must first fit your **null phenotype model** in R/Python. The goal is to extract the inputs your model implies, such as **posterior weights** (for mixture models) and **residuals** (for normal or generalized models).

**Example A: Normal model (single-component)**
```r
fit <- lm(y ~ age + sex + batch, data = pheno)
resid <- resid(fit)
# Single-component model has no mixing: weights are all 1.
weights <- rep(1, length(resid))
```

**Example B: Two-component mixture model (mixtools)**
```r
library(mixtools)

mix <- normalmixEM(y, k = 2)
posterior <- mix$posterior
# Posterior probability of component 1 for each sample
weights <- posterior[, 1]
# Model-implied mean for each sample
mu_hat <- posterior[, 1] * mix$mu[1] + posterior[, 2] * mix$mu[2]
resid <- y - mu_hat
```

You then use these outputs (or other model-derived quantities) to build your **weight matrix** `W` for Stage 2. See `make_weights_*` helpers for common cases, or use `make_weights_generic(...)` for custom models.

**Important:** **Mean-center non-intercept columns of `W` (e.g., multiplicative/variance columns)** to reduce collinearity with the intercept and avoid power loss. This is a re-parameterization (it does not change the model fit), but it can dramatically improve identifiability.

### Step 2: Understanding the “Cutoff”
In a mixture model, there is **no hard cutoff**. The relevant quantity is the **posterior probability** (the weight) computed in Step 1. FungWAS is designed to use these **continuous weights** rather than a hard dichotomy.

If a binary cutoff is required for interpretation, the natural threshold is the point where:

$$P(\text{Class A}) = P(\text{Class B})$$

But for analysis, fungwas prefers the continuous posterior probabilities.

### Step 3: Run Stage 1 & Stage 2
**Stage 1 (Quantile GWAS):**
```bash
fungwas-stage1 \
  --bgen data.bgen \
  --pheno pheno.txt \
  --pheno-col phenotype \
  --covar covar.txt \
  --covar-col age sex batch \
  --out results/chr22
```

**Stage 2 (Parametric mapping):**
```r
W <- make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))

results <- param_gwas_from_file(
  stage1_file = "results/chr22.stage1.tsv.gz",
  W = W,
  se_mode = "tau_cov",
  cov_file = "results/chr22.cov.gz"
)
```

## Core Functions

### `param_gwas_from_file()`
Streamlined wrapper to run Stage 2 directly from Stage 1 output files.

```r
results <- param_gwas_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  W = my_weight_matrix,
  se_mode = "tau_cov",
  cov_file = "path/to/stage1.cov.gz"
)
```

If your covariance file predates `.cov.ids.tsv.gz`, you must opt in to unsafe mode:

```r
results <- param_gwas_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  W = my_weight_matrix,
  se_mode = "tau_cov",
  cov_file = "path/to/stage1.cov.gz",
  unsafe_skip_cov_ids = TRUE
)
```

Covariance storage compatibility:
- Stage 2 auto-detects and reads both formats:
- Legacy float32 `.cov.gz`
- Compact int8 `.cov.gz` with per-SNP scales in `.cov.scale.gz`

### Parameter covariance (generic, any model)
If you run with `se_mode = "tau_cov"`, Stage 2 can emit per-SNP covariance of the
mapped parameters (for mashr or other downstream tools).

```r
out <- param_gwas_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  W = my_weight_matrix,
  se_mode = "tau_cov",
  cov_file = "path/to/stage1.cov.gz",
  return_cov = "upper",
  param_cov_file = "stage2.paramcov.gz"
)

# Upper-triangle covariance per SNP (P x K*(K+1)/2)
param_cov_upper <- out$param_cov_upper
```

For mashr-style inputs:

```r
fit <- param_gwas(
  stage1_object,
  transform = "custom_W",
  transform_args = list(W = my_weight_matrix),
  se_mode = "tau_cov",
  cov_data = my_cov_data,
  return_cov = "mashr"
)

mashr_inputs <- fit$mashr
# mashr_inputs$Bhat  (P x K), mashr_inputs$Shat (P x K), mashr_inputs$V (K x K x P)
```

### `param_gwas()`
Low-level function working on in-memory objects (useful for simulations).

```r
stage2_fit <- param_gwas(stage1_object, transform = "custom_W", transform_args = list(W=W))
```

### `param_gwas_multi()` (model tournament)
Fit multiple models in one pass and return per-SNP AIC/BIC for model selection.

```r
W_list <- list(
  add = matrix(1, nrow = length(taus), ncol = 1),
  add_mult = cbind(1, seq(-1, 1, length.out = length(taus)))
)

model_scores <- param_gwas_multi(stage1_object, W_list, N = 20000)
head(model_scores)
```

### `param_gwas_multi_from_file()`
Convenience wrapper that reads Stage 1 and (optionally) covariance files.

```r
model_scores <- param_gwas_multi_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  cov_file = "path/to/stage1.cov.gz",
  W_list = W_list,
  N = 20000
)
```

If your covariance file predates `.cov.ids.tsv.gz`, you must opt in to unsafe mode:

```r
model_scores <- param_gwas_multi_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  cov_file = "path/to/stage1.cov.gz",
  W_list = W_list,
  N = 20000,
  unsafe_skip_cov_ids = TRUE
)
```

### Weight Builders
Functions to generate the `W` matrix ($K \times T$) mapping $T$ quantiles to $K$ parameters.

-   `make_weights_vqtl(taus, q_tau, mu, sd)`
-   `make_weights_normal_mixture(taus, q_tau, p1, mu1, sd1, ...)`
-   `make_weights_generic(...)`

## Outputs and Definitions

The Stage 2 output includes per-SNP parameter estimates and diagnostics. For score/LRT-style inference, the key statistics are:

- **beta_mu:** Estimated effect of the SNP on the mean of the target distribution parameterization.
- **beta_mu_se:** Standard error of `beta_mu`.
- **beta_sigma2:** Estimated effect of the SNP on the log-variance (variance parameter) of the target distribution.
- **beta_sigma2_se:** Standard error of `beta_sigma2`.
- **Q:** Heterogeneity (misfit) statistic, computed as a quadratic form of residual quantile slopes, testing whether the SNP affects *any* parameters of the chosen model. Under standard large-sample assumptions this is approximately $\chi^2_{df}$ with $df = T - K$.
- **Q_pval:** P-value for the Q statistic under the $\chi^2$ reference distribution.

## Testing
Proof / Validation scripts are located in `tests/`. These are standalone R scripts verifying parameter recovery on simulated data.

```bash
cd tests
Rscript test-vQTL.R
```

## Debugging Weight Matrices
If results look underpowered or unstable, run the identifiability diagnostics on your `W` matrix:

```bash
Rscript stage2/diagnose_identifiability.R
```

Quick wrapper (no examples, optional empirical covariance):

```bash
stage2/diagnose_weights.sh path/to/weights.rds [stage1.tsv.gz] [cov.gz]
```

To explore the effect of different tau grids (optional, advanced):

```bash
Rscript stage2/diagnose_tau_grid.R --weights path/to/weights.rds
```

Notes:
- If you provide `stage1.tsv.gz` and `cov.gz`, the diagnostics use **empirical** tau covariance (recommended).
- If not provided, diagnostics use **theoretical** (normal) tau covariance.

## Generic Tau-Grid Model Comparison
Use this when you already have Stage 1 beta/se matrices and multiple model weight matrices, and you want an analogous model-tournament summary across tau-grid sizes.

Script:

```bash
Rscript stage2/scripts/taugrid_model_compare_fixed_nocov.R \
  --stage1-tsv path/to/stage1.tsv.gz \
  --weights add=path/add.rds,shift=path/shift.rds,both=path/both.rds,vqtl=path/vqtl.rds \
  --out-dir results/taugrid_generic \
  --kmax 40 \
  --k-eval 9,13,17,25,33,40 \
  --greedy-mode both
```

Key options:
- `--greedy-mode=no_cov` uses fixed-monotone greedy selection with no covariance pre-whitening.
- `--greedy-mode=emp_cov` uses fixed-monotone greedy selection after empirical tau-covariance pre-whitening from Stage 1 beta profiles.
- `--greedy-mode=both` runs both and writes direct comparison outputs.

Main outputs:
- `greedy_fixed_monotone_<mode>_sequence.tsv`
- `model_compare_summary_by_k_<mode>.tsv`
- `model_compare_detail_by_snp_<mode>.tsv.gz`
- `winner_counts_by_k_<mode>.tsv`
- `agreement_vs_k_<mode>.png`
- `delta_boxplot_by_k_<mode>.png`

## Known-Truth Simulation for Tau Count and Tau-Covariance
Use this to test whether increasing taus helps/hurts model recovery, and whether covariance-aware scoring (`tau_cov`) improves recovery over diagonal-only (`dwls`) under realistic tau correlation.

Script:

```bash
Rscript stage2/scripts/simulate_taugrid_model_recovery.R \
  --weights add=path/add.rds,shift=path/shift.rds,both=path/both.rds,vqtl=path/vqtl.rds \
  --out-dir results/sim_taugrid \
  --n-snps 400 \
  --n-reps 10 \
  --k-eval 9,13,17,25,33,40 \
  --grid-mode both \
  --rho 0.95 \
  --noise-sd 0.04 \
  --cov-jitter 0.25
```

Design notes:
- Simulates known-truth SNPs from your supplied models (`add`, `shift`, `both`, `vqtl` by default).
- Adds per-SNP correlated tau-noise with AR(1)-like structure and SNP-to-SNP covariance heterogeneity.
- Scores each SNP under all models with both `dwls` and `tau_cov`.
- Reports recovery accuracy by `k`, tau-grid mode, and true model.

Main outputs:
- `sim_summary.tsv`
- `sim_accuracy_by_truth.tsv`
- `sim_confusion_counts.tsv`
- `sim_detail.tsv.gz`
- `sim_accuracy_vs_k.png`
- `sim_accuracy_by_truth_vs_k.png`
- `sim_delta_vs_k.png`

## Full Stage1->Stage2 Recovery Simulation
This workflow runs a true end-to-end benchmark:
1. Simulate raw `G` and `Y`
2. Run Stage1 (with per-SNP tau covariance in `.cov.gz`)
3. Run Stage2 model recovery across tau-grid sizes and `dwls` vs `tau_cov`
4. Evaluate against known truth (`true_model`, `theta`, true tau profile)

### Step 1: Simulate raw data and run Stage1
```bash
python3 stage2/scripts/simulate_full_pipeline_stage1.py \
  --out-prefix results/fullsim/bmi_fullsim \
  --weights add=inputs/data/weights/bmi_100tau/weights_additive_bmi.rds,shift=inputs/data/weights/bmi_100tau/weights_multiplicative_bmi.rds,both=inputs/data/weights/bmi_100tau/weights_both_bmi.rds,vqtl=inputs/data/weights/bmi_100tau/weights_vqtl_bmi.rds \
  --n-samples 12000 \
  --n-snps 1500 \
  --n-threads 8
```

Block/taus stability note:
- Ensure `n_blocks > n_taus` (enforced by Stage 1).
- Recommended: `n_blocks >= n_taus + 20`.
- Practical default tradeoff: 45 interior taus (`0.02-0.98`) with auto blocks `max(n_taus + 20, 64)`.

Outputs:
- `results/fullsim/bmi_fullsim.stage1.tsv.gz`
- `results/fullsim/bmi_fullsim.cov.gz`
- `results/fullsim/bmi_fullsim.cov.ids.tsv.gz`
- `results/fullsim/bmi_fullsim.truth.tsv.gz`
- `results/fullsim/bmi_fullsim.meta.json`

### Step 2: Evaluate Stage2 recovery vs tau count
```bash
Rscript stage2/scripts/evaluate_full_pipeline_taucount.R \
  --stage1-tsv=results/fullsim/bmi_fullsim.stage1.tsv.gz \
  --cov-file=results/fullsim/bmi_fullsim.cov.gz \
  --truth-tsv=results/fullsim/bmi_fullsim.truth.tsv.gz \
  --weights=add=inputs/data/weights/bmi_100tau/weights_additive_bmi.rds,shift=inputs/data/weights/bmi_100tau/weights_multiplicative_bmi.rds,both=inputs/data/weights/bmi_100tau/weights_both_bmi.rds,vqtl=inputs/data/weights/bmi_100tau/weights_vqtl_bmi.rds \
  --out-dir=results/fullsim/bmi_fullsim.eval \
  --k-eval=9,13,17,25,33,40,60 \
  --grid-mode=both
```

### One-command wrapper
```bash
stage2/scripts/run_full_pipeline_simulation.sh \
  results/fullsim/bmi_fullsim \
  add=inputs/data/weights/bmi_100tau/weights_additive_bmi.rds,shift=inputs/data/weights/bmi_100tau/weights_multiplicative_bmi.rds,both=inputs/data/weights/bmi_100tau/weights_both_bmi.rds,vqtl=inputs/data/weights/bmi_100tau/weights_vqtl_bmi.rds \
  12000 1500 8 9,13,17,25,33,40,60
```

Key evaluation outputs:
- `summary_recovery.tsv`: accuracy, delta, profile-MSE by `grid`, `k`, `se_mode`
- `confusion_counts.tsv`: full confusion table (`true_model` vs `winner`)
- `theta_bias_mse.tsv`: bias/MSE/RMSE of true-model parameter estimates
- `accuracy_vs_k.png`: recovery vs tau count
- `winner_mse_vs_k.png`: profile-MSE vs tau count

## REGENIE inputs (fast Stage 1)

If you want to run Stage 1 with REGENIE, use `prepare_regenie()` to generate
RIF phenotypes and covariates.

```r
library(fungwasStage2)
library(data.table)

pheno <- fread("inputs/regenie_pheno_shrink_0.7.txt")
covar <- fread("inputs/regenie_covar.txt")

prepare_regenie(
  pheno_df = pheno,
  covar_df = covar,
  phenotypes = c("testosterone_shrink_0.7_male"),
  taus = seq(0.1, 0.9, 0.05),
  covar_cols = setdiff(names(covar), c("FID", "IID")),
  out_inputs = "regenie_inputs",
  out_rtau = "rtau"
)
```
