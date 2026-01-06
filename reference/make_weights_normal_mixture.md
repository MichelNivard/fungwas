# Build tau-by-K weights for a 2-component Normal mixture

Constructs a weight matrix \\W\\ that maps RIF tau-slopes \\b(\tau)\\
into component-parameter effects under a two-Normal mixture baseline.

## Usage

``` r
make_weights_normal_mixture(
  taus,
  q_tau,
  p1,
  mu1,
  sd1,
  mu2,
  sd2,
  include_membership = FALSE,
  tiny = 1e-12
)
```

## Arguments

- taus:

  Numeric vector of quantile levels.

- q_tau:

  Numeric vector of baseline quantiles at `taus` (`type=8` recommended).

- p1:

  Proportion of class 1 in (0,1); class 2 is `1 - p1`.

- mu1, sd1:

  Mean and standard deviation of component 1.

- mu2, sd2:

  Mean and standard deviation of component 2.

- include_membership:

  Logical; include a first column `"gamma"` if TRUE.

- tiny:

  Small positive floor for stabilizing divisions (default: mixture pdf
  clamped to `tiny`).

## Value

A numeric matrix (`T x K`) with column names:

- If `include_membership = FALSE`: `c("beta_1","beta_2")`

- Else: `c("gamma","beta_1","beta_2")`

## Details

Let \\Y \sim p_1 N(\mu_1, \sigma_1^2) + (1-p_1) N(\mu_2, \sigma_2^2)\\.
For baseline unconditional quantiles \\q\_\tau\\ and mixture pdf
\\f(q\_\tau)\\, the default columns of \\W\\ correspond to component
means: \$\$W_1(\tau) = \frac{p_1 f_1(q\_\tau)}{f(q\_\tau)}, \quad
W_2(\tau) = \frac{(1-p_1) f_2(q\_\tau)}{f(q\_\tau)}\$\$

If `include_membership = TRUE`, a first column is added for membership
(log-odds) perturbation: \$\$W\_\gamma(\tau) = - \frac{p_1
(1-p_1)}{f(q\_\tau)} \left\\F_1(q\_\tau) - F_2(q\_\tau)\right\\\$\$

## Examples

``` r
taus <- seq(0.10, 0.90, by = 0.05)
y <- rnorm(2000, 2, 1)
q_tau <- as.numeric(quantile(y, taus, type = 8))
W <- make_weights_normal_mixture(
  taus, q_tau,
  p1 = 0.5, mu1 = 1.2, sd1 = 0.45,
  mu2 = 3.0, sd2 = 0.7,
  include_membership = TRUE
)
dim(W); colnames(W)
#> [1] 17  3
#> [1] "gamma"  "beta_1" "beta_2"
```
