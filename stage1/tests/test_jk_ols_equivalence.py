"""
Regression test: jackknife OLS vs analytic OLS.

This keeps the original proof logic, but in a small deterministic pytest-sized
case so it can run in CI.
"""

import numpy as np

from fungwas_stage1 import core
from fungwas_stage1.core import (
    compute_block_scores,
    compute_covariance_from_blocks,
    residualize,
)


def simulate_data(n=1200, m=24, seed=42, effect_scale=100.0):
    """Simulate genotype and phenotype data with a small amount of missingness."""
    rng = np.random.default_rng(seed)

    maf = rng.uniform(0.05, 0.5, m)
    g = rng.binomial(2, maf, (n, m)).astype(float)
    g[rng.random((n, m)) < 0.02] = np.nan

    beta_true = rng.normal(size=m) / effect_scale
    g_mean_imputed = np.where(np.isnan(g), np.nanmean(g, axis=0), g)
    y = rng.normal(size=n) + g_mean_imputed @ beta_true
    covariates = rng.normal(size=(n, 2))

    return g, y, covariates


def analytic_ols(y, g, covariates=None, min_mac=5.0, min_var=1e-8):
    """Analytic OLS using the same preprocessing as the Stage 1 kernel."""
    n, m = g.shape

    if covariates is None:
        x = np.ones((n, 1))
    else:
        x = np.column_stack([np.ones(n), covariates])

    y_resid, q = residualize(y, x)
    betas = np.full(m, np.nan)
    ses = np.full(m, np.nan)

    for j in range(m):
        g_j = g[:, j].astype(float).copy()
        valid = ~np.isnan(g_j)
        n_valid = int(np.sum(valid))
        if n_valid < 100:
            continue

        g_sum = np.nansum(g_j)
        g_mean = g_sum / n_valid
        mac = min(g_sum, 2 * n_valid - g_sum)
        if mac < min_mac:
            continue

        g_j[~valid] = g_mean
        g_resid = g_j - q @ (q.T @ g_j)
        g2 = np.sum(g_resid ** 2)
        if g2 < min_var:
            continue

        beta = g_resid @ y_resid / g2
        resid = y_resid - g_resid * beta
        hc0_sum = np.sum((g_resid ** 2) * (resid ** 2))
        se = np.sqrt(hc0_sum) / g2

        betas[j] = beta
        ses[j] = se

    return betas, ses


def jackknife_ols(y, g, covariates=None, n_blocks=32, seed=42, n_threads=1):
    """Jackknife OLS using the stage1 block-score implementation."""
    n, m = g.shape

    if covariates is None:
        x = np.ones((n, 1))
    else:
        x = np.column_stack([np.ones(n), covariates])

    y_resid, q = residualize(y, x)

    rng = np.random.default_rng(seed)
    block_ids = rng.integers(0, n_blocks, n, dtype=np.int32)

    stats = compute_block_scores(
        g,
        y_resid[:, None],
        q,
        block_ids,
        n_blocks,
        n_threads=n_threads,
    )

    betas = np.full(m, np.nan)
    ses = np.full(m, np.nan)
    for j in range(m):
        if np.isnan(stats[j, 0, 0]):
            continue
        beta, se, _ = compute_covariance_from_blocks(stats[j], n_blocks)
        betas[j] = beta[0]
        ses[j] = se[0]

    return betas, ses


def test_jackknife_ols_matches_analytic_ols():
    g, y, covariates = simulate_data()

    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.set_numpy_fallback_allowed(True)
        betas_ols, ses_ols = analytic_ols(y, g, covariates=covariates)
        betas_jk, ses_jk = jackknife_ols(y, g, covariates=covariates, n_blocks=32, n_threads=1)
    finally:
        core.ALLOW_NUMPY_FALLBACK = original_allow

    valid_beta = np.isfinite(betas_ols) & np.isfinite(betas_jk)
    valid_se = np.isfinite(ses_ols) & np.isfinite(ses_jk) & (ses_ols > 0) & (ses_jk > 0)

    assert int(valid_beta.sum()) >= 10
    assert int(valid_se.sum()) >= 10

    beta_corr = np.corrcoef(betas_ols[valid_beta], betas_jk[valid_beta])[0, 1]
    se_ratio = ses_jk[valid_se] / ses_ols[valid_se]

    assert beta_corr > 0.999
    assert np.median(se_ratio) > 0.80
    assert np.median(se_ratio) < 1.20
    assert np.mean(se_ratio) > 0.80
    assert np.mean(se_ratio) < 1.20
