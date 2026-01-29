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

### Step 2: Understanding the “Cutoff” (Alfred’s Question)
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

To explore the effect of different tau grids (optional, advanced):

```bash
Rscript stage2/diagnose_tau_grid.R --weights path/to/weights.rds
```

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
