# Generic weight builder for parametric GWAS

Constructs a tau-by-K weight matrix \\W\\ for an arbitrary parametric
family, mapping RIF tau-slopes \\b(tau)\\ into parameter effects.

## Usage

``` r
make_weights_generic(
  taus,
  q_tau,
  dist_cdf,
  dist_pdf,
  params,
  grad_funcs,
  tiny = 1e-12
)
```

## Arguments

- taus:

  Numeric vector of quantile levels (length T).

- q_tau:

  Numeric vector of baseline quantiles at `taus` (length T).

- dist_cdf:

  Function (y, params) -\> numeric. CDF of the distribution.

- dist_pdf:

  Function (y, params) -\> numeric. PDF of the distribution.

- params:

  Named list of baseline parameters (e.g., `list(mu = 0, sd = 1)`).

- grad_funcs:

  Named list of gradient functions. Each element is a function (y,
  params) -\> numeric giving derivative of F(y; params) wrt that
  parameter. Can be analytic, or constructed with
  [`make_fd_grad()`](https://michelnivard.github.io/fungwas/reference/make_fd_grad.md).

- tiny:

  Small positive floor to stabilize divisions (\\f(q_tau)\\ clamped).

## Value

A T x K numeric matrix of weights with columns named by parameter.

## Details

The generic approach uses the identity: \$\$W_j(tau) = -
(dF/dtheta_j)(q_tau; params) / f(q_tau; params),\$\$ where \\F\\ is the
CDF and \\f\\ is the PDF.

## Examples

``` r
taus <- seq(0.1, 0.9, 0.1)
y <- rnorm(5000, 2, 1)
q_tau <- as.numeric(quantile(y, taus, type = 8))

dist_cdf <- function(y, params) pnorm(y, mean = params$mu, sd = params$sd)
dist_pdf <- function(y, params) dnorm(y, mean = params$mu, sd = params$sd)
grad_mu  <- function(y, params) -dnorm((y - params$mu)/params$sd)/params$sd
grad_sd  <- function(y, params) { z <- (y - params$mu)/params$sd; -(z*dnorm(z))/params$sd }

W <- make_weights_generic(
  taus, q_tau, dist_cdf, dist_pdf,
  params = list(mu = 2, sd = 1),
  grad_funcs = list(beta_mu = grad_mu, beta_sd = grad_sd)
)
```
