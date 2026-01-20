#!/usr/bin/env python3
"""
FungWas Stage 1 CLI

Production-ready RIF quantile GWAS with jackknife covariance estimation.

Usage:
    fungwas-stage1 --bgen data.bgen --pheno pheno.txt --out results/chr22
    
Outputs:
    {out}.stage1.tsv.gz  - Beta/SE per SNP per tau
    {out}.cov.gz         - Upper triangle covariance matrices (binary)
    {out}.Rtau.tsv       - Population tau correlation matrix (for plugin_cor)
"""

import argparse
import gzip
import logging
import os
import struct
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Optional

import numpy as np

from .core import (
    compute_rif,
    residualize,
    compute_block_scores,
    compute_block_scores_multi,
    compute_covariance_from_blocks,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(message)s",
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)


def get_memory_mb() -> float:
    """Get current memory usage in MB."""
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    except:
        return 0.0


def log_mem(msg: str):
    """Log message with memory usage."""
    logger.info(f"[MEM {get_memory_mb():.0f}MB] {msg}")


def read_table_whitespace(file_path: str) -> list[dict]:
    """Read whitespace-delimited tables (space or tab) into list of dicts."""
    rows = []
    with open(file_path) as f:
        header = None
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if header is None:
                header = parts
                continue
            if len(parts) < len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def load_phenotypes(pheno_file: str, pheno_col: str, covar_file: Optional[str], 
                    taus: np.ndarray, keep_file: Optional[str] = None) -> tuple:
    """
    Load raw phenotype and covariates, compute RIF, and residualize.
    
    Parameters
    ----------
    pheno_file : str
        Path to phenotype file (space-delimited with FID, IID, and phenotype column)
    pheno_col : str
        Name of the phenotype column to analyze
    covar_file : str, optional
        Path to covariate file
    taus : array
        Quantile levels for RIF computation
    keep_file : str, optional
        File with sample IDs to keep (one per line, FID or FID\tIID format)
        
    Returns
    -------
    sample_ids : list
        FIDs of included samples
    RIF_resid : array (N, T)
        Residualized RIF matrix
    Q : array (N, K)
        Covariate QR basis  
    q_tau : array (T,)
        Empirical quantiles
    """
    # Load keep list if provided
    keep_set = None
    if keep_file:
        log_mem(f"Loading sample keep list from {keep_file}")
        keep_set = set()
        with open(keep_file) as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    keep_set.add(parts[0])  # Use FID
        log_mem(f"Keeping {len(keep_set)} samples")
    
    log_mem(f"Loading phenotypes from {pheno_file}")
    
    pheno_data = read_table_whitespace(pheno_file)
    
    # Check phenotype column exists
    if pheno_col not in pheno_data[0]:
        available = [c for c in pheno_data[0].keys() if c not in ['FID', 'IID']]
        raise ValueError(f"Phenotype column '{pheno_col}' not found. Available: {available}")
    
    # Load covariates if provided
    if covar_file:
        log_mem(f"Loading covariates from {covar_file}")
        covar_data = read_table_whitespace(covar_file)
        covar_map = {r['FID']: r for r in covar_data}
        covar_cols = [c for c in covar_data[0].keys() if c not in ['FID', 'IID']]
    else:
        covar_map = {}
        covar_cols = []
    
    # Extract phenotype and covariates
    sample_ids = []
    Y_raw = []
    X_list = []
    
    for row in pheno_data:
        fid = row['FID']
        
        # Filter by keep list if provided
        if keep_set is not None and fid not in keep_set:
            continue
        
        try:
            y_val = float(row[pheno_col]) if row[pheno_col] != 'NA' else np.nan
            
            if np.isnan(y_val):
                continue
            
            if covar_file:
                if fid not in covar_map:
                    continue
                c_row = covar_map[fid]
                x_vals = [float(c_row[c]) if c_row[c] != 'NA' else np.nan for c in covar_cols]
                if np.any(np.isnan(x_vals)):
                    continue
                X_list.append(x_vals)
            
            sample_ids.append(fid)
            Y_raw.append(y_val)
            
        except (ValueError, KeyError):
            continue
    
    Y = np.array(Y_raw)
    N = len(Y)
    log_mem(f"Loaded {N} samples with valid phenotype")
    
    # Compute RIF from raw phenotype
    log_mem(f"Computing RIF matrix for {len(taus)} taus...")
    RIF, q_tau, fhat = compute_rif(Y, taus)
    
    # Setup covariate matrix
    if covar_file:
        X = np.column_stack([np.ones(N), np.array(X_list)])
    else:
        X = np.ones((N, 1))
    
    # Residualize RIF on covariates
    log_mem(f"Residualizing RIF on {X.shape[1]} covariates...")
    Q, _ = np.linalg.qr(X)
    RIF_resid = RIF - Q @ (Q.T @ RIF)
    
    return sample_ids, RIF_resid, Q, q_tau


def load_phenotypes_multi(pheno_file: str, pheno_cols: list[str], covar_file: Optional[str],
                          taus: np.ndarray, keep_file: Optional[str] = None) -> tuple:
    """
    Load multiple phenotypes and compute residualized RIF matrices.
    """
    keep_set = None
    if keep_file:
        log_mem(f"Loading sample keep list from {keep_file}")
        keep_set = set()
        with open(keep_file) as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    keep_set.add(parts[0])
        log_mem(f"Keeping {len(keep_set)} samples")

    log_mem(f"Loading phenotypes from {pheno_file}")

    pheno_data = read_table_whitespace(pheno_file)

    for col in pheno_cols:
        if col not in pheno_data[0]:
            available = [c for c in pheno_data[0].keys() if c not in ['FID', 'IID']]
            raise ValueError(f"Phenotype column '{col}' not found. Available: {available}")

    if covar_file:
        log_mem(f"Loading covariates from {covar_file}")
        covar_data = read_table_whitespace(covar_file)
        covar_map = {r['FID']: r for r in covar_data}
        covar_cols = [c for c in covar_data[0].keys() if c not in ['FID', 'IID']]
    else:
        covar_map = {}
        covar_cols = []

    sample_ids = []
    y_values = {col: [] for col in pheno_cols}
    X_list = []

    for row in pheno_data:
        fid = row['FID']
        if keep_set is not None and fid not in keep_set:
            continue

        try:
            pheno_vals = []
            for col in pheno_cols:
                val = float(row[col]) if row[col] != 'NA' else np.nan
                pheno_vals.append(val)
            if np.any(np.isnan(pheno_vals)):
                continue

            if covar_file:
                if fid not in covar_map:
                    continue
                c_row = covar_map[fid]
                x_vals = [float(c_row[c]) if c_row[c] != 'NA' else np.nan for c in covar_cols]
                if np.any(np.isnan(x_vals)):
                    continue
                X_list.append(x_vals)

            sample_ids.append(fid)
            for col, val in zip(pheno_cols, pheno_vals):
                y_values[col].append(val)

        except (ValueError, KeyError):
            continue

    N = len(sample_ids)
    log_mem(f"Loaded {N} samples with valid phenotypes")

    if covar_file:
        X = np.column_stack([np.ones(N), np.array(X_list)])
    else:
        X = np.ones((N, 1))

    log_mem(f"Residualizing phenotypes on {X.shape[1]} covariates...")
    Q, _ = np.linalg.qr(X)

    RIF_resid_list = []
    q_tau_list = []

    for col in pheno_cols:
        log_mem(f"Computing RIF matrix for {col}...")
        t_rif = time.perf_counter()
        Y = np.array(y_values[col])
        RIF, q_tau, _ = compute_rif(Y, taus)
        log_mem(f"RIF build for {col}: {time.perf_counter() - t_rif:.2f}s")
        RIF_resid = RIF - Q @ (Q.T @ RIF)
        RIF_resid_list.append(RIF_resid)
        q_tau_list.append(q_tau)

    return sample_ids, RIF_resid_list, Q, q_tau_list


def extract_snps_bgenix(bgen_path: str, snp_file: str, bgenix_path: str,
                         out_dir: str) -> str:
    """Extract SNPs using bgenix, return path to temp BGEN."""
    temp_bgen = tempfile.NamedTemporaryFile(
        suffix='.bgen', delete=False, dir=out_dir
    ).name
    
    log_mem(f"Extracting SNPs with bgenix...")
    cmd = [bgenix_path, '-g', bgen_path, '-incl-rsids', snp_file]
    
    with open(temp_bgen, 'wb') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
    
    if result.returncode != 0:
        os.unlink(temp_bgen)
        raise RuntimeError(f"bgenix failed: {result.stderr.decode()}")
    
    return temp_bgen


def main():
    parser = argparse.ArgumentParser(
        description="FungWas Stage 1: Fast RIF quantile GWAS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Required arguments
    parser.add_argument('--bgen', required=True,
                        help='Path to BGEN file')
    parser.add_argument('--pheno', required=True,
                        help='Path to phenotype file (space-delimited, with FID IID and phenotype column)')
    parser.add_argument('--pheno-col', dest='pheno_col', required=False,
                        help='Name of the phenotype column to analyze')
    parser.add_argument('--pheno-cols', dest='pheno_cols', default=None,
                        help='Comma-separated phenotype columns to analyze in one pass')
    parser.add_argument('--out', required=True,
                        help='Output prefix (creates {out}.stage1.tsv.gz, {out}.cov.gz)')
    
    # Optional arguments
    parser.add_argument('--sample', default=None,
                        help='Path to sample file (if not embedded in BGEN)')
    parser.add_argument('--covar', default=None,
                        help='Path to covariate file')
    parser.add_argument('--snps', default=None,
                        help='File with SNP IDs to process (one per line)')
    parser.add_argument('--keep', default=None,
                        help='File with sample IDs to keep (one per line)')
    parser.add_argument('--loco-preds', dest='loco_preds', default=None,
                        help='REGENIE step1 LOCO predictions (for future use)')
    
    # Algorithm parameters
    parser.add_argument('--taus', default="0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90",
                        help='Comma-separated tau values')
    parser.add_argument('--blocks', type=int, default=25,
                        help='Number of jackknife blocks')
    parser.add_argument('--batch-size', dest='batch_size', type=int, default=500,
                        help='SNPs per batch for memory efficiency')
    parser.add_argument('--threads', type=int, default=1,
                        help='OpenMP threads for C++ kernel')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for block assignment')
    
    # Tool paths
    parser.add_argument('--bgenix-path', dest='bgenix_path', default='bgenix',
                        help='Path to bgenix binary')
    
    # Output options
    parser.add_argument('--output-rtau', dest='output_rtau', action='store_true',
                        help='Also output R_tau correlation matrix for plugin_cor')
    
    args = parser.parse_args()

    if not args.pheno_col and not args.pheno_cols:
        parser.error("Provide --pheno-col or --pheno-cols")
    
    # Parse taus
    taus = np.array([float(t) for t in args.taus.split(',')])
    T = len(taus)
    
    log_mem(f"Starting FungWas Stage 1 (taus={T}, blocks={args.blocks})")
    t_start = time.time()
    t_load = time.perf_counter()
    
    if args.pheno_cols:
        pheno_cols = [c.strip() for c in args.pheno_cols.split(',') if c.strip()]
    else:
        pheno_cols = [args.pheno_col]

    use_multi = len(pheno_cols) > 1
    if not use_multi and args.pheno_cols:
        log_mem("Single phenotype provided via --pheno-cols; using single-phenotype path.")

    if not use_multi:
        sample_ids, RIF_resid, Q, q_tau = load_phenotypes(
            args.pheno, pheno_cols[0], args.covar, taus, keep_file=args.keep
        )
        RIF_resid_list = [RIF_resid]
        q_tau_list = [q_tau]
    else:
        log_mem(f"Loading {len(pheno_cols)} phenotypes in one pass...")
        sample_ids, RIF_resid_list, Q, q_tau_list = load_phenotypes_multi(
            args.pheno, pheno_cols, args.covar, taus, keep_file=args.keep
        )
    log_mem(f"Phenotype load time: {time.perf_counter() - t_load:.2f}s")
    N = len(sample_ids)
    
    # Setup BGEN
    target_bgen = args.bgen
    temp_bgen = None
    
    if args.snps:
        out_dir = os.path.dirname(args.out) or '.'
        os.makedirs(out_dir, exist_ok=True)
        temp_bgen = extract_snps_bgenix(args.bgen, args.snps, args.bgenix_path, out_dir)
        target_bgen = temp_bgen
    
    # Open BGEN
    from bgen.reader import BgenFile
    bfile = BgenFile(target_bgen)
    bgen_samples = bfile.samples
    
    # Match samples
    sample_to_idx = {s: i for i, s in enumerate(bgen_samples)}
    matched_indices = [sample_to_idx.get(fid) for fid in sample_ids]
    valid_mask = np.array([idx is not None for idx in matched_indices])
    bgen_indices = np.array([idx for idx in matched_indices if idx is not None])
    
    # Filter to matched samples
    RIF_resid_list = [rif[valid_mask] for rif in RIF_resid_list]
    Q = Q[valid_mask]
    N_matched = len(bgen_indices)
    
    log_mem(f"Matched {N_matched}/{N} samples to BGEN")
    t_batches = time.perf_counter()
    
    # Assign jackknife blocks
    np.random.seed(args.seed)
    block_ids = np.random.randint(0, args.blocks, N_matched).astype(np.int32)
    
    # Output setup
    out_dir = os.path.dirname(args.out) or '.'
    os.makedirs(out_dir, exist_ok=True)
    
    def _sanitize_name(name: str) -> str:
        return ''.join(c if c.isalnum() or c in ('-', '_') else '_' for c in name)

    out_prefixes = []
    for col in pheno_cols:
        if len(pheno_cols) == 1:
            out_prefixes.append(args.out)
        else:
            out_prefixes.append(f"{args.out}.{_sanitize_name(col)}")

    out_tsv_list = [f"{prefix}.stage1.tsv.gz" for prefix in out_prefixes]
    out_cov_list = [f"{prefix}.cov.gz" for prefix in out_prefixes]
    
    # Column headers
    cols = ['snp_id', 'chr', 'bp', 'a1', 'a2']
    for t in taus:
        cols.extend([f'beta_tau_{t:.2f}', f'se_tau_{t:.2f}'])
    
    # Upper triangle size
    n_cov_elements = T * (T + 1) // 2
    
    # For R_tau computation
    all_betas_list = [([] if args.output_rtau else None) for _ in pheno_cols]
    
    # Process in batches
    total_snps = 0
    G_batch = np.zeros((N_matched, args.batch_size), dtype=np.float64, order='F')
    batch_info = []
    
    from contextlib import ExitStack

    total_compute = 0.0
    total_write = 0.0

    with ExitStack() as stack:
        f_tsv_list = [stack.enter_context(gzip.open(path, 'wt')) for path in out_tsv_list]
        f_cov_list = [stack.enter_context(gzip.open(path, 'wb')) for path in out_cov_list]
        for f_tsv in f_tsv_list:
            f_tsv.write('\t'.join(cols) + '\n')
        
        idx_in_batch = 0
        
        for variant in bfile:
            # Extract dosages
            probs = variant.probabilities
            dosages = probs[bgen_indices, 1] + 2 * probs[bgen_indices, 2]
            G_batch[:, idx_in_batch] = dosages
            
            batch_info.append((
                variant.rsid, variant.chrom, variant.pos,
                variant.alleles[0], variant.alleles[1]
            ))
            
            idx_in_batch += 1
            
            # Process batch
            if idx_in_batch >= args.batch_size:
                t_batch = time.time()
                
                t_compute = time.perf_counter()
                if use_multi:
                    stats_list = compute_block_scores_multi(
                        G_batch, RIF_resid_list, Q, block_ids, args.blocks,
                        n_threads=args.threads
                    )
                else:
                    stats_list = [
                        compute_block_scores(
                            G_batch, RIF_resid_list[0], Q, block_ids, args.blocks,
                            n_threads=args.threads
                        )
                    ]
                total_compute += time.perf_counter() - t_compute

                t_write = time.perf_counter()
                for p_idx, stats in enumerate(stats_list):
                    f_tsv = f_tsv_list[p_idx]
                    f_cov = f_cov_list[p_idx]
                    all_betas = all_betas_list[p_idx]

                    for i, (snp_stats, info) in enumerate(zip(stats, batch_info)):
                        if np.isnan(snp_stats[0, 0]):
                            continue

                        beta, se, Sigma = compute_covariance_from_blocks(snp_stats, args.blocks)

                        row = list(info)
                        for b, s in zip(beta, se):
                            row.extend([f'{b:.6g}', f'{s:.6g}'])
                        f_tsv.write('\t'.join(map(str, row)) + '\n')

                        cov_upper = Sigma[np.triu_indices(T)]
                        f_cov.write(struct.pack(f'{n_cov_elements}f', *cov_upper))

                        if args.output_rtau:
                            all_betas.append(beta)
                total_write += time.perf_counter() - t_write
                
                total_snps += args.batch_size
                elapsed = time.time() - t_batch
                log_mem(f"Processed {total_snps} SNPs ({args.batch_size/elapsed:.1f} SNPs/s)")
                
                idx_in_batch = 0
                batch_info = []
        
        # Final partial batch
        if idx_in_batch > 0:
            t_compute = time.perf_counter()
            if use_multi:
                stats_list = compute_block_scores_multi(
                    G_batch[:, :idx_in_batch], RIF_resid_list, Q, block_ids, args.blocks,
                    n_threads=args.threads
                )
            else:
                stats_list = [
                    compute_block_scores(
                        G_batch[:, :idx_in_batch], RIF_resid_list[0], Q, block_ids, args.blocks,
                        n_threads=args.threads
                    )
                ]
            total_compute += time.perf_counter() - t_compute

            t_write = time.perf_counter()
            for p_idx, stats in enumerate(stats_list):
                f_tsv = f_tsv_list[p_idx]
                f_cov = f_cov_list[p_idx]
                all_betas = all_betas_list[p_idx]

                for i, (snp_stats, info) in enumerate(zip(stats, batch_info)):
                    if np.isnan(snp_stats[0, 0]):
                        continue

                    beta, se, Sigma = compute_covariance_from_blocks(snp_stats, args.blocks)

                    row = list(info)
                    for b, s in zip(beta, se):
                        row.extend([f'{b:.6g}', f'{s:.6g}'])
                    f_tsv.write('\t'.join(map(str, row)) + '\n')

                    cov_upper = Sigma[np.triu_indices(T)]
                    f_cov.write(struct.pack(f'{n_cov_elements}f', *cov_upper))

                    if args.output_rtau:
                        all_betas.append(beta)
            total_write += time.perf_counter() - t_write
            
            total_snps += idx_in_batch
    
    # Cleanup temp file
    if temp_bgen:
        os.unlink(temp_bgen)
    
    # Output R_tau if requested
    if args.output_rtau:
        for prefix, all_betas in zip(out_prefixes, all_betas_list):
            if not all_betas:
                continue
            log_mem(f"Computing R_tau correlation matrix for {prefix}...")
            beta_mat = np.array(all_betas)
            R_tau = np.corrcoef(beta_mat.T)

            rtau_file = f"{prefix}.Rtau.tsv"
            with open(rtau_file, 'w') as f:
                f.write('\t'.join([f'tau_{t:.2f}' for t in taus]) + '\n')
                for row in R_tau:
                    f.write('\t'.join([f'{v:.6f}' for v in row]) + '\n')
            logger.info(f"Wrote R_tau: {rtau_file}")
    
    elapsed = time.time() - t_start
    log_mem(
        f"Batch compute time: {total_compute:.1f}s, write time: {total_write:.1f}s, "
        f"batch total: {time.perf_counter() - t_batches:.1f}s"
    )
    log_mem(f"Done! {total_snps} SNPs in {elapsed:.1f}s ({total_snps/elapsed:.1f} SNPs/s)")
    logger.info(f"Output: {', '.join(out_tsv_list)}")


if __name__ == "__main__":
    main()
