
# fungwas

Flexible functional GWAS via RIF Quantile GWAS


## Overview

`fungwas` implements a **fast quantile GWAS pipeline** based on **Recentered Influence Functions (RIFs)**.  
It allows you to estimate SNP effects on *distributional parameters* (e.g. means, variances, or mixture components) rather than only on the mean.

Core features:
- RIF-based quantile GWAS with **closed-form OLS slopes** (no `rq()`, no GLS).
- Delta-method standard errors for mapped parameters.
- Generic **tau → parameter mappings** via weight matrices (`W`).
- Built-in weight constructors:
  - `make_weights_normal_mixture()` — for a two-Normal mixture model with SNP effects on both means and mixture membership.
  - `make_weight_vqtl()` — for SNP effects on mean and variance of a Normally distributed phenotype.
- Support for custom `W` matrices, enabling user-defined parametric systems.
- `qgwas_rif()` to run 1. the rif quantile the GWAS and 2. transform the results.

---

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("MichelNivard/fungwas")
```


## Example

```r
library(fungwas)

set.seed(42)
N <- 5000; P <- 50
taus <- seq(0.10, 0.90, by = 0.05)

# Simulate genotype and phenotype
maf <- runif(P, 0.05, 0.5)
G <- matrix(rbinom(N * P, 2, rep(maf, each = N)), N, P)
z  <- rbinom(N, 1, 0.5)
Y  <- rnorm(N, mean = ifelse(z==1, 3.0, 1.2), sd = ifelse(z==1, 0.7, 0.45))

# Run with two_normal transform
fit <- qgwas_rif(
  Y, G, taus = taus,
  transform = "two_normal",
  transform_args = list(
    p1 = 0.5, mu1 = 1.2, sd1 = 0.45,
    mu2 = 3.0, sd2 = 0.7,
    include_membership = TRUE
  ),
  se_mode = "diagonal"
)

# SNP effects on mixture means
head(fit$params)
```


## Functions

* **Main GWAS analysis**

  * `qgwas_rif()` — run RIF quantile GWAS and map the tau-slopes into model parameter estimates.

* **Weight builders**

  * `make_weights_two_normal()` — map tau-slopes into effects on two Normal components (with optional membership).
  * `make_weight_variance()` — map tau-slopes into SNP effects on mean and variance of a Normal.

* **Helpers (internal)**

  * `.build_rif_matrix()` — build RIF(Y; tau).
  * `.residualize_on_C()` — regress SNP or outcome on covariates.
  * `.nearPD_eig()` — repair correlation/covariance matrices.


## License

MIT License © [Your Name]

