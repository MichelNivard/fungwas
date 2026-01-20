# FungWas Stage 2: Parametric Mapping (R Package)

This R package implements "Stage 2" of the FungWas framework: mapping quantile GWAS results (from Stage 1) to parameters of a structural model (like a mixture model or variance-heterogeneity model).

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

## Core Functions

### `param_gwas_from_file()`
Streamlined wrapper to run Stage 2 directly from Stage 1 output files.

```r
results <- param_gwas_from_file(
  stage1_results = "path/to/stage1.tsv.gz",
  W = my_weight_matrix,
  se_mode = "tau_cov",
  cov_file = "path/to/stage1.cov.gz"
)
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

## Testing
Proof / Validation scripts are located in `tests/`. These are standalone R scripts verifying parameter recovery on simulated data.

```bash
cd tests
Rscript test-vQTL.R
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
