"""
Core Stage 1 computation functions.

Provides:
- RIF computation (Recentered Influence Functions)
- Block score accumulation for jackknife
- Covariance estimation from block statistics
"""

import numpy as np
from typing import Tuple, Optional
import logging

logger = logging.getLogger(__name__)

# Try to import C++ extension
try:
    from . import _stage1_cpp as cpp_ext
    HAVE_CPP = True
    logger.info("C++ extension loaded successfully")
except ImportError:
    HAVE_CPP = False
    logger.warning("C++ extension not found, using numpy fallback (slower)")


def resolve_covar_cols(requested: Optional[list[str]], available: list[str]) -> list[str]:
    """
    Resolve requested covariate columns against available columns.

    Parameters
    ----------
    requested : list of str or None
        Covariate columns requested by the user.
    available : list of str
        Available covariate columns in the file (excluding FID/IID).

    Returns
    -------
    list of str
        Columns to use as covariates.
    """
    if requested is None:
        return available
    missing = [c for c in requested if c not in available]
    if missing:
        raise ValueError(f"Covariate column(s) not found: {missing}. Available: {available}")
    return requested


def compute_rif(Y: np.ndarray, taus: np.ndarray, 
                density_floor: float = 1e-8) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute Recentered Influence Function matrix.
    
    Parameters
    ----------
    Y : array (N,)
        Phenotype values
    taus : array (T,)
        Quantile levels
    density_floor : float
        Minimum density value to avoid division issues
        
    Returns
    -------
    RIF : array (N, T)
        RIF values for each sample and tau
    q_tau : array (T,)
        Empirical quantiles at each tau
    fhat : array (T,)
        Kernel density estimates at each quantile
    """
    N = len(Y)
    T = len(taus)
    
    # Empirical quantiles (type 8 for unbiased)
    q_tau = np.quantile(Y, taus, method='median_unbiased')
    
    # Kernel density at each quantile
    # Silverman's rule for bandwidth
    h = 1.06 * np.std(Y) * N ** (-0.2)
    
    fhat = np.zeros(T)
    for t, q in enumerate(q_tau):
        fhat[t] = np.mean(np.exp(-0.5 * ((Y - q) / h) ** 2)) / (h * np.sqrt(2 * np.pi))
    fhat = np.maximum(fhat, density_floor)
    
    # RIF = q + (tau - I(Y <= q)) / f(q)
    RIF = np.zeros((N, T))
    for t in range(T):
        indicator = (Y <= q_tau[t]).astype(float)
        RIF[:, t] = q_tau[t] + (taus[t] - indicator) / fhat[t]
    
    return RIF, q_tau, fhat


def residualize(Y: np.ndarray, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Residualize Y on covariates X using QR decomposition.
    
    Parameters
    ----------
    Y : array (N,) or (N, T)
        Outcome(s) to residualize
    X : array (N, K)
        Covariate matrix (should include intercept)
        
    Returns
    -------
    Y_resid : array, same shape as Y
        Residualized outcome
    Q : array (N, K)
        Orthonormal basis for X
    """
    Q, R = np.linalg.qr(X)
    if Y.ndim == 1:
        Y_resid = Y - Q @ (Q.T @ Y)
    else:
        Y_resid = Y - Q @ (Q.T @ Y)
    return Y_resid, Q


def compute_block_scores(G: np.ndarray, Y_resid: np.ndarray, Q: np.ndarray,
                         block_ids: np.ndarray, n_blocks: int,
                         min_mac: float = 5.0, min_var: float = 1e-8,
                         n_threads: int = 1) -> np.ndarray:
    """
    Compute per-block sufficient statistics for jackknife.
    
    Parameters
    ----------
    G : array (N, M)
        Genotype dosages (may contain NaN for missing)
    Y_resid : array (N, T)
        Residualized RIF matrix
    Q : array (N, K)
        QR basis for covariates
    block_ids : array (N,)
        Block assignment for each sample (0 to n_blocks-1)
    n_blocks : int
        Number of jackknife blocks
    min_mac : float
        Minimum minor allele count filter
    min_var : float
        Minimum variance filter
    n_threads : int
        Number of OpenMP threads (if C++ available)
        
    Returns
    -------
    stats : array (M, n_blocks, 1+T)
        Per-SNP, per-block statistics
        [:, :, 0] = sum(G_resid^2) in block
        [:, :, 1:] = sum(G_resid * Y_resid) in block per tau
    """
    if HAVE_CPP:
        return cpp_ext.stage1_block_scores(
            np.asfortranarray(G, dtype=np.float64),
            np.asfortranarray(Y_resid, dtype=np.float64),
            np.asfortranarray(Q, dtype=np.float64),
            block_ids.astype(np.int32),
            n_blocks, min_mac, min_var, n_threads
        )
    else:
        return _compute_block_scores_numpy(G, Y_resid, Q, block_ids, n_blocks, min_mac, min_var)


def compute_block_scores_multi(G: np.ndarray, RIF_resid_list: list[np.ndarray], Q: np.ndarray,
                               block_ids: np.ndarray, n_blocks: int,
                               min_mac: float = 5.0, min_var: float = 1e-8,
                               n_threads: int = 1) -> list[np.ndarray]:
    """
    Compute per-block sufficient statistics for multiple phenotypes.

    Parameters
    ----------
    G : array (N, M)
        Genotype dosages (may contain NaN for missing)
    RIF_resid_list : list of arrays (N, T)
        Residualized RIF matrices
    Q : array (N, K)
        QR basis for covariates
    block_ids : array (N,)
        Block assignment for each sample (0 to n_blocks-1)
    n_blocks : int
        Number of jackknife blocks
    min_mac : float
        Minimum minor allele count filter
    min_var : float
        Minimum variance filter
    n_threads : int
        Number of OpenMP threads (if C++ available)

    Returns
    -------
    list of arrays (M, n_blocks, 1+T)
        Per-SNP, per-block statistics for each phenotype
    """
    if not isinstance(RIF_resid_list, (list, tuple)):
        raise ValueError("RIF_resid_list must be a list of arrays")

    if HAVE_CPP and hasattr(cpp_ext, "stage1_block_scores_multi"):
        return cpp_ext.stage1_block_scores_multi(
            np.asfortranarray(G, dtype=np.float64),
            [np.asfortranarray(rif, dtype=np.float64) for rif in RIF_resid_list],
            np.asfortranarray(Q, dtype=np.float64),
            block_ids.astype(np.int32),
            n_blocks, min_mac, min_var, n_threads
        )

    return [
        _compute_block_scores_numpy(G, rif, Q, block_ids, n_blocks, min_mac, min_var)
        for rif in RIF_resid_list
    ]


def _compute_block_scores_numpy(G, Y_resid, Q, block_ids, n_blocks, min_mac, min_var):
    """Pure numpy fallback (slower)."""
    N, M = G.shape
    T = Y_resid.shape[1]
    
    stats = np.full((M, n_blocks, 1 + T), np.nan)
    
    for j in range(M):
        g = G[:, j].copy()
        valid = ~np.isnan(g)
        n_valid = np.sum(valid)
        
        if n_valid < 100:
            continue
            
        g_mean = np.nanmean(g)
        mac = min(np.nansum(g), 2 * n_valid - np.nansum(g))
        
        if mac < min_mac:
            continue
            
        # Mean impute missing
        g[~valid] = g_mean
        
        # Residualize on covariates
        g_resid = g - Q @ (Q.T @ g)
        
        g2 = np.sum(g_resid ** 2)
        if g2 < min_var:
            continue
        
        # Accumulate per block
        block_stats = np.zeros((n_blocks, 1 + T))
        for i in range(N):
            b = block_ids[i]
            if b >= n_blocks:
                continue
            block_stats[b, 0] += g_resid[i] ** 2
            block_stats[b, 1:] += g_resid[i] * Y_resid[i, :]
        
        stats[j] = block_stats
    
    return stats


def compute_covariance_from_blocks(stats: np.ndarray, n_blocks: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute jackknife estimates from block statistics.
    
    Parameters
    ----------
    stats : array (n_blocks, 1+T)
        Block statistics for one SNP
    n_blocks : int
        Number of jackknife blocks
        
    Returns
    -------
    beta : array (T,)
        Point estimates
    se : array (T,)  
        Standard errors
    Sigma : array (T, T)
        Full covariance matrix
    """
    T = stats.shape[1] - 1
    D_k = stats[:, 0]  # sum(G^2) per block
    N_k = stats[:, 1:]  # sum(G*Y) per block per tau
    
    D_total = np.sum(D_k)
    N_total = np.sum(N_k, axis=0)
    
    if D_total <= 0:
        return np.full(T, np.nan), np.full(T, np.nan), np.full((T, T), np.nan)
    
    # Leave-one-out betas
    jk_betas = np.zeros((n_blocks, T))
    for k in range(n_blocks):
        denom = D_total - D_k[k]
        if denom > 0:
            jk_betas[k] = (N_total - N_k[k]) / denom
    
    # Point estimate
    beta = N_total / D_total
    
    # Jackknife covariance
    centered = jk_betas - np.mean(jk_betas, axis=0)
    Sigma = ((n_blocks - 1) / n_blocks) * (centered.T @ centered)
    
    se = np.sqrt(np.diag(Sigma))
    
    return beta, se, Sigma


def run_stage1(G: np.ndarray, Y: np.ndarray, 
               taus: np.ndarray,
               covariates: Optional[np.ndarray] = None,
               n_blocks: int = 25,
               seed: int = 42,
               n_threads: int = 1) -> dict:
    """
    Run full Stage 1 quantile GWAS pipeline.
    
    Parameters
    ----------
    G : array (N, M)
        Genotype dosage matrix
    Y : array (N,)
        Phenotype
    taus : array (T,)
        Quantile levels to analyze
    covariates : array (N, K), optional
        Covariate matrix (intercept added automatically)
    n_blocks : int
        Number of jackknife blocks
    seed : int
        Random seed for block assignment
    n_threads : int
        OpenMP threads for C++ kernel
        
    Returns
    -------
    dict with keys:
        betas : array (M, T) - tau-level effect estimates
        se : array (M, T) - standard errors
        covariances : list of (T, T) arrays - full covariance per SNP
        q_tau : array (T,) - empirical quantiles
        taus : array (T,) - tau levels used
    """
    N, M = G.shape
    T = len(taus)
    
    # Compute RIF
    logger.info(f"Computing RIF matrix for {T} taus...")
    RIF, q_tau, fhat = compute_rif(Y, taus)
    
    # Setup covariates
    if covariates is None:
        X = np.ones((N, 1))
    else:
        X = np.column_stack([np.ones(N), covariates])
    
    # Residualize RIF
    logger.info("Residualizing RIF on covariates...")
    RIF_resid, Q = residualize(RIF, X)
    
    # Assign blocks
    np.random.seed(seed)
    block_ids = np.random.randint(0, n_blocks, N).astype(np.int32)
    
    # Compute block scores
    logger.info(f"Computing block scores for {M} SNPs...")
    all_stats = compute_block_scores(G, RIF_resid, Q, block_ids, n_blocks, n_threads=n_threads)
    
    # Aggregate results
    betas = np.zeros((M, T))
    ses = np.zeros((M, T))
    covariances = []
    
    for j in range(M):
        if np.isnan(all_stats[j, 0, 0]):
            betas[j] = np.nan
            ses[j] = np.nan
            covariances.append(np.full((T, T), np.nan))
        else:
            beta, se, Sigma = compute_covariance_from_blocks(all_stats[j], n_blocks)
            betas[j] = beta
            ses[j] = se
            covariances.append(Sigma)
    
    return {
        'betas': betas,
        'se': ses,
        'covariances': covariances,
        'q_tau': q_tau,
        'taus': taus,
        'n_blocks': n_blocks,
    }
