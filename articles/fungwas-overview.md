# Functional GWAS via Quantile Regression and RIF

## Motivation

Most GWAS focus on **mean effects**: does a SNP increase or decrease the
*average* level of a trait?

But genetic effects can be more subtle:

- Some SNPs may increase **variability** (vQTLs).
- Others may shift the balance between **subtypes or mixture
  components**.
- Some may affect **both mean and variance** together.

The `fungwas` package provides a **fast, flexible framework** to test
SNP effects on *distributional parameters*, not just means.

------------------------------------------------------------------------

## Core idea

1.  Use **Recentered Influence Functions (RIFs)** to link SNPs to
    changes in trait quantiles.

    - This gives you “SNP slopes” across the distribution.
    - It is fast (closed-form OLS, no iterative quantile regression).

2.  Combine quantile slopes using a **weight matrix W** to map them into
    any parameter of interest:

    - mean, variance, mixture membership, etc.

3.  Estimate SNP effects and standard errors on those parameters.

This two-step process is the backbone of `fungwas`.

------------------------------------------------------------------------

## Example 1: Standard mean effect

Let’s simulate a simple trait with SNP effects on the mean.

``` r
N <- 5000
P <- 50
maf <- runif(P, 0.05, 0.5)
G <- matrix(rbinom(N * P, 2, rep(maf, each = N)), N, P)
colnames(G) <- paste0("SNP", seq_len(P))

# Simulate a mean effect
beta_true <- rnorm(P) / 40
Y <- rnorm(N, mean = 2 + G %*% beta_true, sd = 1)

taus <- seq(0.1, 0.9, 0.05)

# Stage 1: quantile GWAS
stage1 <- quantile_gwas(Y, G, taus = taus)
#> Building RIF matrix on raw Y...
#> Computing per-SNP tau-slopes...

# Stage 2: mean/variance mapping (vQTL weights)
W <- make_weights_vqtl(taus, stage1$q_tau, mu = mean(Y), sd = sd(Y))
fit <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W))

head(t(fit$params))
#>          beta_mu  beta_sigma2
#> SNP1 -0.03798818 -0.026820633
#> SNP2 -0.01832948  0.018785274
#> SNP3  0.02361505 -0.045563528
#> SNP4  0.04302755  0.001305224
#> SNP5  0.05758468 -0.048737974
#> SNP6  0.01513990  0.199157417
```

Here, the column `beta_mu` corresponds to the SNP effect on the mean.

------------------------------------------------------------------------

## Example 2: Variance effects (vQTLs)

Suppose SNPs affect not just the mean but also the **spread** of the
phenotype.

``` r
beta_mu_true     <- rnorm(P) / 40
beta_sigma2_true <- rnorm(P) / 40

mu <- 2 + G %*% beta_mu_true
sigma2 <- 1 + G %*% beta_sigma2_true
sigma2[sigma2 <= 0] <- 0.1

Y <- rnorm(N, mean = mu, sd = sqrt(sigma2))

taus <- seq(0.05, 0.95, 0.05)
stage1 <- quantile_gwas(Y, G, taus = taus)
#> Building RIF matrix on raw Y...
#> Computing per-SNP tau-slopes...

W_var <- make_weights_vqtl(taus, stage1$q_tau, mu = mean(Y), sd = sd(Y))
fit_var <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W_var))

head(t(fit_var$params))
#>          beta_mu beta_sigma2
#> SNP1 -0.01474599 -0.10471233
#> SNP2  0.03817296  0.02671172
#> SNP3 -0.06506831 -0.05781526
#> SNP4  0.01184463  0.03804201
#> SNP5 -0.01047337 -0.04669658
#> SNP6  0.01771871 -0.08678875
```

Now we obtain SNP effects on both the **mean** and the **variance**.

------------------------------------------------------------------------

## Example 3: Mixture GWAS

Many complex traits are mixtures of underlying subtypes. For example,
“cases” might consist of two symptom clusters, or biomarker
distributions may show multimodality.

We can model SNP effects on **component means** and **class membership**
in a two-Normal mixture.

``` r
p1 <- 0.5
mu1 <- 1.2; sd1 <- 0.5
mu2 <- 3.0; sd2 <- 0.8

# Simulate phenotype
z <- rbinom(N, 1, p1)
Y <- ifelse(z == 1, rnorm(N, mu1, sd1), rnorm(N, mu2, sd2))

taus <- seq(0.1, 0.9, 0.05)
stage1 <- quantile_gwas(Y, G, taus = taus)
#> Building RIF matrix on raw Y...
#> Computing per-SNP tau-slopes...

fit_mix <- param_gwas(
  stage1,
  transform = "two_normal",
  transform_args = list(
    p1 = p1, mu1 = mu1, sd1 = sd1,
    mu2 = mu2, sd2 = sd2,
    include_membership = TRUE
  )
)

head(t(fit_mix$params))
#>            gamma        beta_1        beta_2
#> SNP1 -0.05727378 -1.903029e-02 -0.0233213774
#> SNP2  0.07954904 -7.991301e-05  0.0076335131
#> SNP3 -0.06057390 -2.434824e-02  0.0003425944
#> SNP4  0.02074416 -7.745207e-04 -0.0188630721
#> SNP5 -0.12819602 -1.526194e-02 -0.0227835128
#> SNP6  0.11976032  8.322496e-02 -0.0342396586
```

## Key functions and workflow

The workflow in `fungwas` always follows **two stages**:

1.  **Quantile GWAS**: estimate SNP effects on *quantile slopes* across
    the phenotype distribution.
    - Function:
      [`quantile_gwas()`](https://michelnivard.github.io/fungwas/reference/quantile_gwas.md)  
    - Inputs:
      - `Y`: vector of phenotypes (length N).  
      - `G`: genotype matrix, N × P (rows = individuals, cols = SNPs).  
      - `C`: optional covariates (N × K, e.g. age, sex, PCs).  
      - `taus`: grid of quantile levels (e.g. `seq(0.1, 0.9, 0.05)`).  
    - Output:
      - A list containing RIF slopes per SNP × τ, their SEs, and
        baseline quantiles.  
    - Think of this as a **quantile-level GWAS**.
2.  **Parameter GWAS**: map quantile slopes into parameter effects using
    a weight matrix `W`.
    - Function:
      [`param_gwas()`](https://michelnivard.github.io/fungwas/reference/param_gwas.md)  
    - Inputs:
      - The output of
        [`quantile_gwas()`](https://michelnivard.github.io/fungwas/reference/quantile_gwas.md).  
      - A mapping (`W`) that tells the software how to combine τ-slopes
        into parameter effects.  
      - Either supply a custom `W` or use a pre-built weight
        constructor.  
    - Output:
      - SNP effects on parameters (means, variances, mixture membership,
        etc).  
    - This is the **interpretation stage**, translating distributional
      shifts into biologically meaningful parameters.

### Weight builders

Weight matrices (`W`) define how quantile slopes correspond to parameter
changes.  
`fungwas` provides several ready-made constructors:

- **Variance GWAS (vQTLs)**
  - [`make_weights_vqtl()`](https://michelnivard.github.io/fungwas/reference/make_weights_vqtl.md)
    — use when you want SNP effects on the *mean* and *variance* of a
    Normal trait.  
  - Example: height variability, BMI dispersion.
- **Mixture GWAS**
  - [`make_weights_normal_mixture()`](https://michelnivard.github.io/fungwas/reference/make_weights_normal_mixture.md)
    — use for a two-component Normal mixture with SNP effects on:
    - Component means.  
    - (Optionally) class membership probability.  
  - Example: SNPs shifting balance between subtypes of cases/controls.
- **Mixture vQTL GWAS**
  - `make_weights_mixture_vqtl()` — extended version where SNPs can also
    affect the **component variances** (in addition to means and
    membership).  
  - Example: genetic effects on both subtype prevalence and
    within-subtype variability.
- **Generic system**
  - [`make_weights_generic()`](https://michelnivard.github.io/fungwas/reference/make_weights_generic.md)
    — the most flexible constructor.
    - You provide the distribution’s CDF, PDF, and derivatives wrt
      parameters.  
    - Or use finite-difference helpers
      ([`make_fd_grad()`](https://michelnivard.github.io/fungwas/reference/make_fd_grad.md)).  
  - Example: log-normal GWAS on `meanlog` and `sdlog`.  
  - Use this when your phenotype is better described by a non-standard
    distribution.

### Putting it together: a typical workflow

1.  **Prepare inputs**:

    - Phenotype vector `Y`.  
    - Genotype dosage matrix `G` (SNPs already QC’d).  
    - Optional covariate matrix `C` (e.g. sex, age, ancestry PCs).

2.  **Run quantile GWAS**:

    ``` r
    taus <- seq(0.1, 0.9, 0.05)
    stage1 <- quantile_gwas(Y, G, taus = taus, C = covariates)
    ```

3.  **Choose a weight system**:

    - If mean & variance:

      ``` r
      W <- make_weights_vqtl(taus, stage1$q_tau, mu = mean(Y), sd = sd(Y))
      ```

    - If two-component mixture:

      ``` r
      W <- make_weights_normal_mixture(taus, stage1$q_tau,
                                       p1 = 0.5, mu1 = 1, sd1 = 1,
                                       mu2 = 3, sd2 = 1.5,
                                       include_membership = TRUE)
      ```

    - If custom:

      ``` r
      W <- make_weights_generic(taus, stage1$q_tau, dist_cdf, dist_pdf, params, grad_funcs)
      ```

4.  **Map to parameters**:

    ``` r
    fit <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W))
    ```

5.  **Inspect results**:

    ``` r
    head(t(fit$params))   # SNP effects on chosen parameters
    head(t(fit$SE_params)) # Standard errors
    ```

In practice:

- **Stage 1 (`quantile_gwas`) is run only once** per dataset.
- You can then re-use it with different weight systems (`param_gwas`) to
  test multiple hypotheses (means, variances, mixtures, etc) *without
  rerunning the GWAS*.

## Why use fungwas?

- Go beyond mean effects: test **variance, mixtures, heterogeneity**.
- Fast: closed-form OLS, no heavy quantile regression.
- Flexible: any distributional system can be defined via
  `make_weights_generic`.
- Useful for **vQTL discovery**, **subtype genetics**, and **causal
  inference** when mean effects don’t tell the full story.

## Further reading

- Koenker & Bassett (1978), Quantile Regression.
- Firpo, Fortin, & Lemieux (2009), Influence Function Regression.
- Recent applications of vQTL and mixture genetics in GWAS.

------------------------------------------------------------------------

## Summary

`fungwas` allows you to re-think GWAS: instead of only asking *“does
this SNP affect the mean?”*, you can ask *“does this SNP affect the
distribution?”* — variance, mixtures, or other user-defined systems.
