# Your custom functional GWAS with \`make_weights_generic\`

## Overview

In `fungwas`, a core idea is to map **quantile GWAS slopes** into
**parameters of a parametric system** (means, variances, mixture
components, etc). This mapping is accomplished by **weight matrices**
`W`.

While some specific constructors are built-in
(e.g. `make_weights_normal_mixture`, `make_weights_vqtl`), the function
[`make_weights_generic()`](https://michelnivard.github.io/fungwas/reference/make_weights_generic.md)
allows you to define your own mapping for **any distribution** so long
as you can supply:

- The **CDF** and **PDF** of the distribution.
- The **gradients of the CDF wrt parameters** (analytic or
  finite-difference).

This vignette illustrates its use with a **log-normal phenotype**.

## Example: Log-Normal GWAS

We simulate phenotypes from a log-normal distribution with SNP effects
on both the `meanlog` and `sdlog` parameters.

``` r
# Simulate genotypes
N <- 5000
P <- 50
maf <- runif(P, 0.05, 0.5)
G <- matrix(rbinom(N * P, 2, rep(maf, each = N)), nrow = N, ncol = P)
colnames(G) <- paste0("SNP", seq_len(P))

# True SNP effects
beta_meanlog_true <- rnorm(P) / 40
beta_sdlog_true   <- rnorm(P) / 40

# Baseline parameters
meanlog0 <- 0.5
sdlog0   <- 0.7

# Individual parameters
meanlog <- meanlog0 + G %*% beta_meanlog_true
sdlog   <- pmax(sdlog0 + G %*% beta_sdlog_true, 0.1)

# Phenotype
Y <- rlnorm(N, meanlog = meanlog, sdlog = sdlog)
```

## Step 1: Define distribution and gradients

We define the **CDF** and **PDF** of the log-normal, and supply
gradients of the CDF wrt `meanlog` and `sdlog`. The derivatives
determine how SNP effects on parameters are expressed in quantile space.

``` r
dist_cdf <- function(y, params) plnorm(y, meanlog = params$meanlog, sdlog = params$sdlog)
dist_pdf <- function(y, params) dlnorm(y, meanlog = params$meanlog, sdlog = params$sdlog)

grad_meanlog <- function(y, params) {
  z <- (log(y) - params$meanlog) / params$sdlog
  -dnorm(z) / params$sdlog
}
grad_sdlog <- function(y, params) {
  z <- (log(y) - params$meanlog) / params$sdlog
  -(z * dnorm(z)) / params$sdlog
}
```

## Step 2: Build weights

We compute weights for a grid of quantiles.

``` r
taus <- seq(0.1, 0.9, 0.05)
q_tau <- as.numeric(quantile(Y, taus, type = 8))

W <- make_weights_generic(
  taus, q_tau,
  dist_cdf, dist_pdf,
  params = list(meanlog = meanlog0, sdlog = sdlog0),
  grad_funcs = list(beta_meanlog = grad_meanlog, beta_sdlog = grad_sdlog)
)

head(W)
#>         beta_meanlog beta_sdlog
#> tau0.1     0.8193355 -0.8184712
#> tau0.15    0.9722950 -0.7335216
#> tau0.2     1.1053469 -0.6313754
#> tau0.25    1.2319771 -0.5128188
#> tau0.3     1.3678373 -0.3649568
#> tau0.35    1.4923085 -0.2124957
```

Each column of `W` corresponds to the mapping from quantile slopes into
parameter-space effects.

## Step 3: Run two-stage GWAS

First, run a **quantile GWAS** to estimate SNP effects on each quantile
(RIF slopes). Then, map these into parameter-space using
[`param_gwas()`](https://michelnivard.github.io/fungwas/reference/param_gwas.md)
with our custom `W`.

``` r
stage1 <- quantile_gwas(Y, G, taus = taus, benchmark = FALSE)
#> Building RIF matrix on raw Y...
#> Computing per-SNP tau-slopes...

fit <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W),
  se_mode = "diagonal"
)

head(t(fit$params))
#>      beta_meanlog  beta_sdlog
#> SNP1 -0.069766804  0.03844083
#> SNP2  0.017892104 -0.01579062
#> SNP3 -0.009028385  0.02589287
#> SNP4  0.021612584 -0.02085744
#> SNP5  0.045506929 -0.01453348
#> SNP6  0.013347512 -0.03759752
```

## Step 4: Compare to true effects

We can check recovery by correlating estimated effects with truth.

``` r
est <- t(fit$params)

cor_meanlog <- cor(beta_meanlog_true, est[, "beta_meanlog"])
cor_sdlog   <- cor(beta_sdlog_true,   est[, "beta_sdlog"])

cor_meanlog
#> [1] 0.8096836
cor_sdlog
#> [1] 0.8384027
```

------------------------------------------------------------------------

## Finite-Difference Gradients

Instead of hand-coding gradients, you can generate them automatically
using
[`make_fd_grad()`](https://michelnivard.github.io/fungwas/reference/make_fd_grad.md).

``` r
grad_fd_meanlog <- make_fd_grad("meanlog", dist_cdf)
grad_fd_sdlog   <- make_fd_grad("sdlog", dist_cdf)

W_fd <- make_weights_generic(
  taus, q_tau,
  dist_cdf, dist_pdf,
  params = list(meanlog = meanlog0, sdlog = sdlog0),
  grad_funcs = list(beta_meanlog = grad_fd_meanlog, beta_sdlog = grad_fd_sdlog)
)

# Run GWAS with FD weights
fit_fd <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W_fd))

cor(est[, "beta_meanlog"], t(fit_fd$params)[, "beta_meanlog"])
#> [1] 1
cor(est[, "beta_sdlog"],   t(fit_fd$params)[, "beta_sdlog"])
#> [1] 1
```

Analytic and FD versions give nearly identical results, making FD a
convenient option when analytic gradients are complicated.

------------------------------------------------------------------------

## Summary

- [`make_weights_generic()`](https://michelnivard.github.io/fungwas/reference/make_weights_generic.md)
  enables flexible GWAS of **arbitrary parametric systems**.
- You provide a CDF, PDF, and parameter gradients (either analytic or
  via
  [`make_fd_grad()`](https://michelnivard.github.io/fungwas/reference/make_fd_grad.md)).
- This generalizes RIF-based GWAS to models far beyond simple mean
  effects.
