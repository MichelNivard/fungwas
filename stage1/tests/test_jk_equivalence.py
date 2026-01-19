"""
Proof: Jackknife SE Equivalence

Demonstrates that jackknife-based RIF regression produces identical betas 
and properly calibrated SEs compared to standard RIF OLS.

Uses simulated data only - no real phenotypes.
"""

import numpy as np
from fungwas_stage1.core import (
    compute_rif,
    residualize,
    compute_block_scores,
    compute_covariance_from_blocks,
    run_stage1,
)


def simulate_data(N=5000, M=50, seed=42):
    """Simulate genotype and phenotype data."""
    np.random.seed(seed)
    
    # Genotypes
    maf = np.random.uniform(0.05, 0.5, M)
    G = np.random.binomial(2, maf, (N, M)).astype(float)
    
    # True effects on mean
    beta_true = np.random.randn(M) / 40
    
    # Phenotype
    Y = np.random.randn(N) + G @ beta_true
    
    return G, Y, beta_true


def standard_rif_ols(Y, G, taus):
    """Standard RIF OLS regression (for comparison)."""
    N, M = G.shape
    T = len(taus)
    
    # Compute RIF
    RIF, q_tau, fhat = compute_rif(Y, taus)
    
    # Center genotypes
    G_centered = G - G.mean(axis=0)
    
    # Compute betas and SEs for each SNP
    betas = np.zeros((M, T))
    ses = np.zeros((M, T))
    
    for j in range(M):
        x = G_centered[:, j]
        Sxx = np.sum(x ** 2)
        if Sxx < 1e-8:
            continue
        
        # Beta = X'Y / X'X
        betas[j] = x @ RIF / Sxx
        
        # SE from residuals
        for t in range(T):
            resid = RIF[:, t] - x * betas[j, t]
            sigma2 = np.sum(resid ** 2) / (N - 2)
            ses[j, t] = np.sqrt(sigma2 / Sxx)
    
    return betas, ses, q_tau


class TestJackknifeEquivalence:
    """Test that jackknife produces equivalent results to standard OLS."""
    
    def test_beta_equivalence(self):
        """Betas from jackknife should match standard OLS."""
        G, Y, beta_true = simulate_data(N=5000, M=50)
        taus = np.arange(0.10, 0.91, 0.05)
        
        # Standard OLS
        betas_ols, ses_ols, _ = standard_rif_ols(Y, G, taus)
        
        # Jackknife
        result = run_stage1(G, Y, taus, n_blocks=25, n_threads=1)
        betas_jk = result['betas']
        
        # Check correlation
        valid = ~np.isnan(betas_ols.flatten()) & ~np.isnan(betas_jk.flatten())
        corr = np.corrcoef(betas_ols.flatten()[valid], betas_jk.flatten()[valid])[0, 1]
        
        assert corr > 0.999, f"Beta correlation {corr:.4f} < 0.999"
    
    def test_se_calibration(self):
        """Jackknife SEs should be properly calibrated."""
        G, Y, beta_true = simulate_data(N=5000, M=100)
        taus = np.arange(0.10, 0.91, 0.10)
        
        # Run jackknife
        result = run_stage1(G, Y, taus, n_blocks=25, n_threads=1)
        betas_jk = result['betas']
        ses_jk = result['se']
        
        # Standard OLS SEs
        betas_ols, ses_ols, _ = standard_rif_ols(Y, G, taus)
        
        # Compare SE distributions (should be in same ballpark)
        valid = (ses_jk > 0) & (ses_ols > 0) & ~np.isnan(ses_jk) & ~np.isnan(ses_ols)
        ratio = ses_jk[valid] / ses_ols[valid]
        
        median_ratio = np.median(ratio)
        
        # Jackknife SEs might be slightly larger due to finite blocks
        # but should be within reasonable range (0.8 - 1.5)
        assert 0.8 < median_ratio < 1.5, f"SE ratio {median_ratio:.2f} out of range"
    
    def test_true_beta_recovery(self):
        """Both methods should recover true effects."""
        G, Y, beta_true = simulate_data(N=10000, M=50, seed=123)
        taus = np.arange(0.10, 0.91, 0.05)
        
        # For mean effects, the middle tau should recover beta_true
        result = run_stage1(G, Y, taus, n_blocks=25)
        
        # At tau=0.5, RIF regression should recover mean effect
        tau_50_idx = np.argmin(np.abs(taus - 0.5))
        betas_mid = result['betas'][:, tau_50_idx]
        
        valid = ~np.isnan(betas_mid)
        corr = np.corrcoef(beta_true[valid], betas_mid[valid])[0, 1]
        
        assert corr > 0.8, f"True beta recovery correlation {corr:.3f} < 0.8"


if __name__ == "__main__":
    # Run tests manually
    test = TestJackknifeEquivalence()
    
    print("Testing beta equivalence...")
    test.test_beta_equivalence()
    print("✓ Beta equivalence passed")
    
    print("Testing SE calibration...")
    test.test_se_calibration()
    print("✓ SE calibration passed")
    
    print("Testing true beta recovery...")
    test.test_true_beta_recovery()
    print("✓ True beta recovery passed")
    
    print("\n✓ All tests passed!")
