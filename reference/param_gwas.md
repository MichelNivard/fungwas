# Parametric GWAS mapping

Maps quantile-level GWAS results into effects on user-specified
parametric system parameters via a weight matrix `W`.

## Usage

``` r
param_gwas(
  stage1,
  transform = c("custom_W", "two_normal"),
  transform_args = list(),
  se_mode = c("diagonal", "plugin_cor", "dwls"),
  plugin_R = NULL
)
```

## Arguments

- stage1:

  Output list from
  [`quantile_gwas()`](https://michelnivard.github.io/fungwas/reference/quantile_gwas.md).

- transform:

  One of:

  - `"custom_W"` – supply a weight matrix `W` in `transform_args`.

  - `"two_normal"` – construct weights for a two-component normal
    mixture.

- transform_args:

  Arguments for the chosen transform.

- se_mode:

  Standard error mode:

  - `"diagonal"` – assume independence across taus (fast).

  - `"plugin_cor"` – plugin correlation with near-PD repair.

  - `"dwls"` – diagonal weighted least squares (Q test calibrated).

- plugin_R:

  Optional T x T correlation matrix for `"plugin_cor"`.

## Value

A list with the following elements:

- `W`:

  Weight matrix used in the parametric mapping.

- `A`:

  Mapping matrix (tau-by-parameter).

- `params`:

  Estimated parameter values.

- `SE_params`:

  Delta-method standard errors of parameter estimates.

- `Q`:

  Vector of Q statistics per SNP (if computed).

- `df`:

  Degrees of freedom for Q statistic.
