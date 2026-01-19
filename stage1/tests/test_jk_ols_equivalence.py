"""
Proof 1: Jackknife SE Equivalence (Standard OLS)

This test validates the block-jackknife SEs for ordinary least squares
regression (Y ~ G) using the exact same preprocessing as the C++ kernel:

- mean imputation of missing genotypes
- MAC/min-variance filtering
- residualization on covariates via QR

This is the first check in the sequential validation chain:
1) OLS jackknife vs analytic OLS
2) RIF/quantile jackknife vs analytic RIF OLS (stage1)
3) Stage2 SE modes vs RIF regression (stage2)
"""

import numpy as np
from fungwas_stage1.core import (
    residualize,
    compute_block_scores,
    compute_covariance_from_blocks,
)


def simulate_data(N=5000, M=50, seed=42, effect_scale=100):
    """Simulate genotype and phenotype data."""
    np.random.seed(seed)

    maf = np.random.uniform(0.05, 0.5, M)
    G = np.random.binomial(2, maf, (N, M)).astype(float)

    beta_true = np.random.randn(M) / effect_scale
    Y = np.random.randn(N) + G @ beta_true

    return G, Y, beta_true


def analytic_ols(Y, G, covariates=None, min_mac=5.0, min_var=1e-8):
    """Analytic OLS using the same preprocessing as the C++ kernel."""
    N, M = G.shape

    if covariates is None:
        X = np.ones((N, 1))
    else:
        X = np.column_stack([np.ones(N), covariates])

    Y_resid, Q = residualize(Y, X)
    betas = np.full(M, np.nan)
    ses = np.full(M, np.nan)

    for j in range(M):
        g = G[:, j].astype(float).copy()
        valid = ~np.isnan(g)
        n_valid = np.sum(valid)
        if n_valid < 100:
            continue

        g_sum = np.nansum(g)
        g_mean = g_sum / n_valid
        mac = min(g_sum, 2 * n_valid - g_sum)
        if mac < min_mac:
            continue

        g[~valid] = g_mean
        g_resid = g - Q @ (Q.T @ g)
        g2 = np.sum(g_resid ** 2)
        if g2 < min_var:
            continue

        beta = g_resid @ Y_resid / g2
        resid = Y_resid - g_resid * beta
        hc0_sum = np.sum((g_resid ** 2) * (resid ** 2))
        se = np.sqrt(hc0_sum) / g2

        betas[j] = beta
        ses[j] = se

    return betas, ses, Y_resid, Q


def jackknife_ols(Y, G, covariates=None, n_blocks=25, seed=42, n_threads=1):
    """Jackknife OLS using the stage1 block-score implementation."""
    N, M = G.shape

    if covariates is None:
        X = np.ones((N, 1))
    else:
        X = np.column_stack([np.ones(N), covariates])

    Y_resid, Q = residualize(Y, X)

    np.random.seed(seed)
    block_ids = np.random.randint(0, n_blocks, N).astype(np.int32)

    stats = compute_block_scores(
        G,
        Y_resid[:, None],
        Q,
        block_ids,
        n_blocks,
        n_threads=n_threads,
    )

    betas = np.full(M, np.nan)
    ses = np.full(M, np.nan)
    for j in range(M):
        if np.isnan(stats[j, 0, 0]):
            continue
        beta, se, _ = compute_covariance_from_blocks(stats[j], n_blocks)
        betas[j] = beta[0]
        ses[j] = se[0]

    return betas, ses


def run_large_scale_proof():
    print("Running OLS jackknife proof (N=100,000)...")
    import matplotlib.pyplot as plt
    G, Y, _ = simulate_data(N=100000, M=100, effect_scale=100)

    print("Running analytic OLS...")
    betas_ols, ses_ols, _, _ = analytic_ols(Y, G)

    print("Running jackknife OLS...")
    betas_jk, ses_jk = jackknife_ols(Y, G, n_blocks=25, n_threads=16)

    valid = np.isfinite(betas_ols) & np.isfinite(betas_jk)
    beta_corr = np.corrcoef(betas_ols[valid], betas_jk[valid])[0, 1]
    print(f"Beta Correlation: {beta_corr:.6f}")

    valid_se = np.isfinite(ses_ols) & np.isfinite(ses_jk) & (ses_ols > 0) & (ses_jk > 0)
    se_ratio = ses_jk[valid_se] / ses_ols[valid_se]

    print("\n--- SE Comparison Stats ---")
    print(f"Median Ratio (JK/OLS): {np.median(se_ratio):.4f}")
    print(f"Mean Ratio   (JK/OLS): {np.mean(se_ratio):.4f}")
    print(f"Std Ratio            : {np.std(se_ratio):.4f}")

    if np.any(valid_se):
        print(f"Analytic SE min      : {np.min(ses_ols[valid_se]):.6f}")
        print(f"Jackknife SE min     : {np.min(ses_jk[valid_se]):.6f}")
        quantiles = np.quantile(ses_ols[valid_se], [0.01, 0.05, 0.5, 0.95, 0.99])
        print("Analytic SE quantiles:", ", ".join(f"{q:.6f}" for q in quantiles))

    slope, intercept = np.polyfit(ses_ols[valid_se], ses_jk[valid_se], 1)
    print(f"Regression slope (SE_JK ~ SE_OLS): {slope:.4f}")

    try:
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        plt.scatter(betas_ols[valid], betas_jk[valid], alpha=0.5, s=8)
        plt.title(f"Betas: JK vs OLS (r={beta_corr:.4f})")
        plt.xlabel("Analytic OLS Beta")
        plt.ylabel("Jackknife Beta")
        plt.plot([betas_ols[valid].min(), betas_ols[valid].max()],
                 [betas_ols[valid].min(), betas_ols[valid].max()],
                 'r--')

        plt.subplot(1, 2, 2)
        plt.scatter(ses_ols[valid_se], ses_jk[valid_se], alpha=0.5, s=8)
        plt.title(f"SEs: JK vs OLS (Slope={slope:.3f})")
        plt.xlabel("Analytic OLS SE")
        plt.ylabel("Jackknife SE")

        mx = max(ses_ols[valid_se].max(), ses_jk[valid_se].max())
        mn = min(ses_ols[valid_se].min(), ses_jk[valid_se].min())
        plt.plot([mn, mx], [mn, mx], 'r--', label='1:1')

        x_vals = np.array([mn, mx])
        plt.plot(x_vals, intercept + slope * x_vals, 'g-', label=f'Fit (k={slope:.2f})')
        plt.legend()

        plt.tight_layout()
        outfile = "proof_stage1_ols_jk_vs_ols.png"
        plt.savefig(outfile)
        print(f"\nSaved comparison plot to {outfile}")
    except Exception as e:
        print(f"Could not save plot: {e}")


if __name__ == "__main__":
    run_large_scale_proof()
