# Build tau-by-2 weights for Normal mean/variance perturbations

Constructs a weight matrix W that maps RIF tau-slopes b(tau) into
effects on the mean and variance of a Normal baseline phenotype.

## Usage

``` r
make_weights_vqtl(taus, q_tau, mu, sd)
```

## Arguments

- taus:

  Numeric vector of quantile levels (length T).

- q_tau:

  Numeric vector of baseline quantiles at taus (length T, type=8
  recommended).

- mu:

  Baseline mean of Y.

- sd:

  Baseline standard deviation of Y.

## Value

A T x 2 numeric matrix W with columns:

- `"beta_mu"` effect of SNP on the mean

- `"beta_sigma2"` effect of SNP on the variance

## Details

Baseline: Y ~ N(mu, sigma^2). For tau-quantiles q_tau = mu + sigma \*
z_tau, z_tau = Phi^-1(tau): \$\$ W_mu(tau) = 1, \quad W_sigma2(tau) =
z_tau / (2 \* sigma). \$\$

## Examples

``` r
taus <- seq(0.1, 0.9, by = 0.2)
y <- rnorm(2000, mean = 2, sd = 1.5)
q_tau <- as.numeric(quantile(y, taus, type = 8))
W <- make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))
W
#>        beta_mu beta_sigma2
#> tau0.1       1  -0.4126440
#> tau0.3       1  -0.1688506
#> tau0.5       1   0.0000000
#> tau0.7       1   0.1688506
#> tau0.9       1   0.4126440
```
