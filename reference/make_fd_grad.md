# Construct a finite-difference gradient function

Builds a function that computes the numerical gradient of a
distribution's CDF with respect to a single parameter using central
differences.

## Usage

``` r
make_fd_grad(param, dist_cdf, eps = 1e-06)
```

## Arguments

- param:

  Name of the parameter (string).

- dist_cdf:

  Function (y, params) -\> numeric. CDF of the distribution.

- eps:

  Finite-difference step size (default 1e-6).

## Value

A function with signature (y, params) that returns dF(y; params) /
d(param) for the given parameter.

## Examples

``` r
dist_cdf <- function(y, params) pgamma(y, shape = params$shape, scale = params$scale)
grad_shape <- make_fd_grad("shape", dist_cdf)
grad_shape(1.5, list(shape = 2, scale = 1))
#> [1] -0.3134886
```
