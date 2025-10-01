# fungwas: Flexible functional GWAS via RIF Quantile GWAS

## Overview

`fungwas` implements a **fast quantile GWAS pipeline** based on **Recentered Influence Functions (RIFs)**.  
The results of the quantile GWAS, which might be of interest in and of themselves, enable the second stage estimation of SNP effects on *distributional parameters* (e.g. means, variances, or mixture components) of parametric models of the phenotype, rather than only on the mean.

Core features:

- RIF-based **quantile GWAS** with closed-form OLS slopes.
- Delta-method based standard errors for estimated parameters.
- Generic **tau → parameter mappings** via weight matrices (`W`).
- Built-in weight constructors:

  - `make_weights_normal_mixture()` —  
    Construct weights for a **two-component Normal mixture** with SNP effects on the two component means.  
    Optionally, include an effect on **class membership** (log-odds of being in component 1).  

  - `make_weight_vqtl()` —  
    Construct weights for SNP effects on the **mean** and **variance** of a single Normal distribution (variance QTL analysis).  

  - `make_weights_normal_mixture_vqtl()` —  
    Construct weights for a **two-component Normal mixture with vQTLs**.  
    Supports SNP effects on:
      * mixture membership (γ, log-odds of component 1),  
      * both component means (μ₁, μ₂), and  
      * both component standard deviations (σ₁, σ₂).  
    This enables testing for **mixture vQTLs**, where genetic variants can affect variability differently across mixture components.  

  - `make_weights_generic()` —  
    A **generic gradient-based builder** for arbitrary parametric models.  
    Given a function that returns quantiles for model parameters, this constructs the Jacobian of quantiles with respect to parameters (analytically or via finite differences) and produces the appropriate weight matrix.  
    This allows you to study custom parametric systems such as log-Normal models, skew-Normals, or other user-defined distributions.

- Support for fully custom `W` matrices, enabling user-defined parametric systems.  

- Two-stage workflow:
  1. `quantile_gwas()` — run the RIF-based quantile GWAS (tau-slopes).
  2. `param_gwas()` — map tau-slopes into parameters of a chosen parametric model.

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("MichelNivard/fungwas")
````

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
Y  <- rnorm(N, mean = ifelse(z == 1, 3.0, 1.2),
               sd   = ifelse(z == 1, 0.7, 0.45))

# Stage 1: tau-level GWAS
stage1 <- quantile_gwas(Y, G, taus = taus)

# Stage 2: map tau slopes into mixture parameters
fit <- param_gwas(
  stage1,
  transform = "two_normal",
  transform_args = list(
    p1 = 0.5, mu1 = 1.2, sd1 = 0.45,
    mu2 = 3.0, sd2 = 0.7,
    include_membership = TRUE
  ),
  se_mode = "diagonal"
)

# SNP effects on mixture membership and means
head(fit$params)
```

## Functions

**Main workflow**

* `quantile_gwas()` — run RIF quantile GWAS (tau-slopes).
* `param_gwas()` — map per SNP tau-slopes into user-specified parametric model parameters.

**Weight builders**

These functions construct weight matrices (`W`) that map quantile‐level SNP effects (tau‐slopes) into parameters of specific parametric systems. They provide ready‐to‐use setups for common use cases:

* `make_weights_normal_mixture()` —
  Map tau‐slopes to SNP effects on the means of two Normal components, with optional effects on **class membership** (log‐odds of belonging to component 1).

* `make_weight_vqtl()` —
  Map tau‐slopes to SNP effects on the **mean** and **variance** of a single Normal distribution (variance QTL analysis).

* `make_weights_normal_mixture_vqtl()` —
  Map tau‐slopes to SNP effects on a **two‐component Normal mixture** with full flexibility:

  * class membership (γ),
  * component means (μ₁, μ₂),
  * component standard deviations (σ₁, σ₂).
    This enables mixture vQTL analysis, where genetic variants influence not just means but also variability across mixture components.

* `make_weights_generic()` —
  Build weights for **any user-defined parametric system**. Requires a function that maps model parameters to quantiles; gradients are computed analytically or numerically (via finite differences).

**Internal helpers**

* `.build_rif_matrix()` — build RIF(Y; tau).
* `.residualize_on_C()` — regress SNPs or outcome on covariates.
* `.nearPD_eig()` — repair correlation/covariance matrices.

## License

MIT License © Michel Nivard & Fergus Hamilton
