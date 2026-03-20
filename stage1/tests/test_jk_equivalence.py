"""
Regression test: Stage 1 jackknife RIF path vs analytic RIF OLS.

This is a pytest-sized version of the original proof script.
"""

import numpy as np

from fungwas_stage1 import core
from fungwas_stage1.core import (
    compute_block_scores,
    compute_covariance_from_blocks,
    compute_rif,
    residualize,
    run_stage1,
)


def simulate_data(n=900, m=18, seed=42, effect_scale=100.0):
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


def analytic_rif_ols(y, g, taus, covariates=None, min_mac=5.0, min_var=1e-8):
    """Analytic RIF OLS using the same preprocessing as the Stage 1 kernel."""
    n, m = g.shape
    t = len(taus)

    rif, _, _ = compute_rif(y, taus)

    if covariates is None:
        x = np.ones((n, 1))
    else:
        x = np.column_stack([np.ones(n), covariates])

    rif_resid, q = residualize(rif, x)
    betas = np.full((m, t), np.nan)
    ses = np.full((m, t), np.nan)

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

        betas[j] = g_resid @ rif_resid / g2

        for tau_idx in range(t):
            resid = rif_resid[:, tau_idx] - g_resid * betas[j, tau_idx]
            hc0_sum = np.sum((g_resid ** 2) * (resid ** 2))
            ses[j, tau_idx] = np.sqrt(hc0_sum) / g2

    return betas, ses


def jackknife_rif(y, g, taus, covariates=None, n_blocks=24, seed=42, n_threads=1):
    """Jackknife RIF estimates using the stage1 block-score implementation."""
    n, m = g.shape

    rif, _, _ = compute_rif(y, taus)

    if covariates is None:
        x = np.ones((n, 1))
    else:
        x = np.column_stack([np.ones(n), covariates])

    rif_resid, q = residualize(rif, x)

    np.random.seed(seed)
    block_ids = np.random.randint(0, n_blocks, n).astype(np.int32)

    stats = compute_block_scores(
        g,
        rif_resid,
        q,
        block_ids,
        n_blocks,
        n_threads=n_threads,
    )

    t = len(taus)
    betas = np.full((m, t), np.nan)
    ses = np.full((m, t), np.nan)
    for j in range(m):
        if np.isnan(stats[j, 0, 0]):
            continue
        beta, se, _ = compute_covariance_from_blocks(stats[j], n_blocks)
        betas[j] = beta
        ses[j] = se

    return betas, ses


def test_stage1_jackknife_matches_analytic_rif_ols():
    g, y, covariates = simulate_data()
    taus = np.array([0.10, 0.30, 0.50, 0.70, 0.90])

    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.set_numpy_fallback_allowed(True)
        betas_ols, ses_ols = analytic_rif_ols(y, g, taus, covariates=covariates)
        betas_jk, ses_jk = jackknife_rif(y, g, taus, covariates=covariates, n_blocks=24, n_threads=1)
    finally:
        core.ALLOW_NUMPY_FALLBACK = original_allow

    valid_beta = np.isfinite(betas_ols) & np.isfinite(betas_jk)
    valid_se = np.isfinite(ses_ols) & np.isfinite(ses_jk) & (ses_ols > 0) & (ses_jk > 0)

    assert int(valid_beta.sum()) >= 30
    assert int(valid_se.sum()) >= 30

    beta_corr = np.corrcoef(betas_ols[valid_beta], betas_jk[valid_beta])[0, 1]
    se_ratio = ses_jk[valid_se] / ses_ols[valid_se]

    assert beta_corr > 0.999
    assert np.median(se_ratio) > 0.75
    assert np.median(se_ratio) < 1.25
    assert np.mean(se_ratio) > 0.75
    assert np.mean(se_ratio) < 1.25


def test_run_stage1_matches_jackknife_rif_reference():
    g, y, covariates = simulate_data(seed=7)
    taus = np.array([0.10, 0.30, 0.50, 0.70, 0.90])

    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.set_numpy_fallback_allowed(True)
        ref_betas, ref_ses = jackknife_rif(y, g, taus, covariates=covariates, n_blocks=24, seed=7, n_threads=1)
        result = run_stage1(
            g,
            y,
            taus,
            covariates=covariates,
            n_blocks=24,
            seed=7,
            n_threads=1,
            allow_numpy_fallback=True,
        )
    finally:
        core.ALLOW_NUMPY_FALLBACK = original_allow

    np.testing.assert_allclose(result["betas"], ref_betas, rtol=1e-10, atol=1e-10, equal_nan=True)
    np.testing.assert_allclose(result["se"], ref_ses, rtol=1e-10, atol=1e-10, equal_nan=True)
