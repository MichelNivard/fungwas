# Quantile-level RIF GWAS

Runs a genome-wide association scan at user-specified quantile levels
(taus) using Recentered Influence Functions (RIF). Supports two
covariate modes:

- "residualize": partial out C from SNPs (but never from Y).

- "include": include C directly in the regression alongside each SNP.

## Usage

``` r
quantile_gwas(
  Y,
  G,
  taus = seq(0.1, 0.9, 0.05),
  C = NULL,
  covar_mode = c("residualize", "include"),
  density_floor = 1e-08,
  benchmark = TRUE,
  verbose = TRUE
)
```

## Arguments

- Y:

  Numeric vector (length N), phenotype.

- G:

  Numeric matrix (N x P) of SNP dosages or genotypes.

- taus:

  Numeric vector of quantile levels.

- C:

  Optional N x K covariate matrix.

- covar_mode:

  Either "residualize" (default) or "include".

  - "residualize": partial out C from G before regression.

  - "include": include C directly in regression.

- density_floor:

  Positive scalar; lower bound for estimated densities.

- benchmark:

  Logical; if TRUE, include timing info.

- verbose:

  Logical; print progress messages.

## Value

A list with taus, q_tau, fhat_tau, Q_slope, SE_tau, timing.
