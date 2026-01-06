# Build tau-by-K weights for a two-component Normal mixture with mean, variance, and membership effects

Constructs a weight matrix `W` that maps RIF tau-slopes \\b(tau)\\ into
SNP effects on the parameters of a two-component Normal mixture: \$\$Y ~
p1 \* N(mu1, sd1^2) + (1 - p1) \* N(mu2, sd2^2).\$\$

## Usage

``` r
make_weights_normal_mixture_vqtl(
  taus,
  q_tau,
  p1,
  mu1,
  sd1,
  mu2,
  sd2,
  tiny = 1e-12
)
```

## Arguments

- taus:

  Numeric vector of quantile levels (length T).

- q_tau:

  Numeric vector of baseline quantiles of Y at `taus` (length T,
  `type=8` recommended).

- p1:

  Proportion of component 1 (in (0,1)); component 2 proportion = 1 - p1.

- mu1, sd1:

  Mean and sd of component 1.

- mu2, sd2:

  Mean and sd of component 2.

- tiny:

  Small positive floor to stabilize divisions (mixture pdf clamped to at
  least `tiny`).

## Value

A T x 5 numeric matrix `W` with columns:
`c("gamma","beta_mu1","beta_mu2","beta_sigma1","beta_sigma2")`.

## Details

Supported SNP effect targets:

- `gamma`: log-odds of membership in component 1 (class proportion).

- `beta_mu1`, `beta_mu2`: effects on the component means.

- `beta_sigma1`, `beta_sigma2`: effects on the component standard
  deviations.

For baseline quantiles \\q_tau\\ and mixture pdf \$\$f(q) = p1 f1(q) +
(1 - p1) f2(q),\$\$ where \\f1, f2\\ are the component pdfs and \\F1,
F2\\ their CDFs, the weights are:

- Membership: \$\$W_gamma(tau) = - (p1 \* (1 - p1) / f(q_tau)) \*
  (F1(q_tau) - F2(q_tau))\$\$

- Means: \$\$W_mu1(tau) = (p1 \* f1(q_tau)) / f(q_tau), \quad W_mu2(tau)
  = ((1 - p1) \* f2(q_tau)) / f(q_tau)\$\$

- Standard deviations: \$\$W_sigma1(tau) = (p1 \* f1(q_tau) \* (q_tau -
  mu1) / sd1) / f(q_tau),\$\$ \$\$W_sigma2(tau) = ((1 - p1) \* f2(q_tau)
  \* (q_tau - mu2) / sd2) / f(q_tau).\$\$

These arise from the derivative of the mixture CDF with respect to each
parameter, mapped into the quantile effect scale via the implicit
function theorem.

## Examples

``` r
taus <- seq(0.1, 0.9, 0.1)
set.seed(1)
y <- c(rnorm(500, 1.2, 0.6), rnorm(500, 3.0, 0.9))
q_tau <- as.numeric(quantile(y, taus, type = 8))
W <- make_weights_normal_mixture_vqtl(
  taus, q_tau,
  p1 = 0.5, mu1 = 1.2, sd1 = 0.6,
  mu2 = 3.0, sd2 = 0.9
)
head(W)
#>             gamma  beta_mu1   beta_mu2 beta_sigma1 beta_sigma2
#> tau0.1 -0.2084154 0.9635124 0.03648760  -0.7484892 -0.09187173
#> tau0.2 -0.2684789 0.9432446 0.05675541  -0.3074682 -0.12584448
#> tau0.3 -0.3559549 0.9040579 0.09594215   0.1073999 -0.18428582
#> tau0.4 -0.4982207 0.8155360 0.18446399   0.4986371 -0.29373769
#> tau0.5 -0.7001093 0.6182523 0.38174773   0.7180740 -0.46790617
#> tau0.6 -0.8293443 0.3026462 0.69735385   0.5369821 -0.56983570
```
