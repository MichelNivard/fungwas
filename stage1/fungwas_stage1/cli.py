#!/usr/bin/env python3
"""
FungWas Stage 1 CLI

Production-ready RIF quantile GWAS with jackknife covariance estimation.

Usage:
    fungwas-stage1 --bgen data.bgen --pheno pheno.txt --out results/chr22
    
Outputs:
    {out}.stage1.tsv.gz  - Beta/SE per SNP per tau
    {out}.cov.gz         - Upper triangle covariance matrices (binary)
    {out}.cov.ids.tsv.gz - SNP IDs for covariance alignment (one per row)
    {out}.Rtau.tsv       - Population tau correlation matrix (for plugin_cor)
"""

import argparse
import gzip
import json
import logging
import os
import re
import struct
import sys
import time
from pathlib import Path
from typing import Optional, List

import numpy as np

from .core import (
    compute_rif,
    residualize,
    compute_block_scores,
    compute_block_scores_multi,
    compute_covariance_from_blocks,
    resolve_covar_cols,
    set_numpy_fallback_allowed,
    require_cpp_extension,
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


def parse_list_arg(values: Optional[list[str]]) -> Optional[list[str]]:
    if values is None:
        return None
    parts = []
    for v in values:
        parts.extend([p.strip() for p in v.split(',') if p.strip()])
    return parts or None


def normalize_chr_label(chrom) -> str:
    """Normalize chromosome labels to plain autosome-style strings."""
    chrom_str = str(chrom).strip()
    chrom_str = re.sub(r'^(chr|CHR)', '', chrom_str)
    return chrom_str


def read_keep_ids(keep_file: str) -> set[str]:
    """Load a sample keep list using the first whitespace-delimited column."""
    keep_set = set()
    with open(keep_file) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                keep_set.add(parts[0])
    return keep_set


def build_sample_aliases(fid: str, iid: Optional[str]) -> list[str]:
    """Build plausible sample identifiers seen across phenotype/BGEN/REGENIE files."""
    aliases = []
    if fid:
        aliases.append(fid)
    if iid and iid != fid:
        aliases.append(iid)
    if fid and iid:
        aliases.extend([
            f"{fid}_{iid}",
            f"{fid} {iid}",
            f"{fid}:{iid}",
        ])
    return aliases


def resolve_sample_column(sample_lookup: dict[str, int], aliases: list[str]) -> Optional[int]:
    """Resolve the first matching sample alias in a LOCO header lookup."""
    for alias in aliases:
        if alias in sample_lookup:
            return sample_lookup[alias]
    return None


def resolve_loco_files(loco_path: str, pheno_cols: list[str]) -> dict[str, str]:
    """
    Resolve REGENIE LOCO inputs to a phenotype->.loco path mapping.

    Supports either a direct `.loco` file or a REGENIE `_pred.list` master file.
    """
    path = Path(loco_path)
    if path.suffix == '.loco':
        if len(pheno_cols) != 1:
            raise ValueError(
                "--loco-preds points to a single .loco file but multiple phenotypes were requested. "
                "Provide the REGENIE *_pred.list file instead."
            )
        return {pheno_cols[0]: str(path)}

    entries: list[list[str]] = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                entries.append(parts)

    loco_files = [parts[-1] for parts in entries if parts[-1].endswith('.loco')]
    if not loco_files:
        raise ValueError(f"No .loco files found in {loco_path}")

    if len(pheno_cols) == 1 and len(loco_files) == 1:
        return {pheno_cols[0]: loco_files[0]}

    mapping = {}
    remaining = set(pheno_cols)

    for parts in entries:
        loco_file = parts[-1]
        if not loco_file.endswith('.loco'):
            continue
        joined = ' '.join(parts[:-1] + [Path(loco_file).stem])
        for pheno in list(remaining):
            if pheno in joined:
                mapping[pheno] = loco_file
                remaining.remove(pheno)
                break

    if remaining:
        unresolved = sorted(remaining)
        unassigned = [f for f in loco_files if f not in mapping.values()]
        if len(unassigned) == len(unresolved):
            for pheno, loco_file in zip(unresolved, unassigned):
                mapping[pheno] = loco_file
        else:
            raise ValueError(
                f"Could not map LOCO files in {loco_path} to phenotypes {unresolved}. "
                "Use a pred.list whose entries include phenotype names or pass a single .loco for one phenotype."
            )

    return mapping


def load_loco_predictions_for_file(loco_file: str,
                                   sample_aliases: list[list[str]]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """
    Load per-chromosome REGENIE LOCO predictions aligned to the phenotype sample order.
    """
    with open(loco_file) as f:
        header = f.readline().strip().split()
        if not header or len(header) < 2:
            raise ValueError(f"Malformed LOCO header in {loco_file}")
        sample_lookup = {sample_id: idx for idx, sample_id in enumerate(header[1:])}

        matched_cols = []
        for idx, aliases in enumerate(sample_aliases):
            col_idx = resolve_sample_column(sample_lookup, aliases)
            matched_cols.append(col_idx)
        valid_mask = np.array([idx is not None for idx in matched_cols], dtype=bool)
        missing = np.where(~valid_mask)[0]
        matched_cols_arr = np.array(
            [idx for idx in matched_cols if idx is not None],
            dtype=np.int64
        )

        if matched_cols_arr.size == 0:
            raise ValueError(f"No phenotype samples were found in LOCO file {loco_file}")
        if missing.size > 0:
            preview = ', '.join(sample_aliases[i][0] if sample_aliases[i] else f"sample_{i}" for i in missing[:5])
            logger.warning(
                "%d phenotype samples are absent from LOCO file %s and will be dropped. Examples: %s",
                int(missing.size), loco_file, preview
            )

        chr_to_pred = {}
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            chrom = normalize_chr_label(parts[0])
            values = np.asarray(parts[1:], dtype=np.float64)
            if values.shape[0] != len(header) - 1:
                raise ValueError(
                    f"Chromosome {parts[0]} in {loco_file} has {values.shape[0]} predictions, "
                    f"expected {len(header) - 1}"
                )
            chr_to_pred[chrom] = values[matched_cols_arr]

    if not chr_to_pred:
        raise ValueError(f"No chromosome predictions found in {loco_file}")
    return chr_to_pred, valid_mask


def load_loco_predictions(loco_path: str,
                          pheno_cols: list[str],
                          sample_aliases: list[list[str]]) -> tuple[dict[str, dict[str, np.ndarray]], np.ndarray]:
    """Load aligned LOCO predictions for one or more phenotypes and return the shared valid sample mask."""
    loco_files = resolve_loco_files(loco_path, pheno_cols)
    raw_predictions = {}
    valid_masks = {}
    combined_valid_mask = np.ones(len(sample_aliases), dtype=bool)
    for pheno in pheno_cols:
        loco_file = loco_files[pheno]
        log_mem(f"Loading LOCO predictions for {pheno} from {loco_file}")
        pheno_preds, valid_mask = load_loco_predictions_for_file(loco_file, sample_aliases)
        raw_predictions[pheno] = pheno_preds
        valid_masks[pheno] = valid_mask
        combined_valid_mask &= valid_mask

    if not np.any(combined_valid_mask):
        raise ValueError("No samples remain after intersecting phenotype rows with LOCO predictions.")

    predictions = {}
    for pheno in pheno_cols:
        pheno_valid_idx = np.where(valid_masks[pheno])[0]
        keep_positions = np.nonzero(combined_valid_mask[pheno_valid_idx])[0]
        predictions[pheno] = {
            chrom: values[keep_positions]
            for chrom, values in raw_predictions[pheno].items()
        }

    return predictions, combined_valid_mask


def load_snps_list(snp_file: str) -> List[str]:
    """Load list of SNP rsIDs from file (one per line)."""
    snps = []
    with open(snp_file) as f:
        for line in f:
            rsid = line.strip().split()[0] if line.strip() else ''
            if rsid:
                snps.append(rsid)
    return snps


def load_analysis_inputs(pheno_file: str, pheno_cols: list[str], covar_file: Optional[str],
                         keep_file: Optional[str] = None,
                         covar_cols: Optional[list[str]] = None) -> tuple:
    """Load phenotypes/covariates and return aligned raw inputs for analysis."""
    keep_set = None
    if keep_file:
        log_mem(f"Loading sample keep list from {keep_file}")
        keep_set = read_keep_ids(keep_file)
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
        available_covar_cols = [c for c in covar_data[0].keys() if c not in ['FID', 'IID']]
        covar_cols = resolve_covar_cols(covar_cols, available_covar_cols)
    elif covar_cols is not None:
        raise ValueError("covar_cols provided but no covariate file (--covar) was supplied.")
    else:
        covar_map = {}
        covar_cols = []

    sample_ids = []
    sample_aliases = []
    y_values = {col: [] for col in pheno_cols}
    X_list = []

    for row in pheno_data:
        fid = row['FID']
        iid = row.get('IID', fid)
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
            sample_aliases.append(build_sample_aliases(fid, iid))
            for col, val in zip(pheno_cols, pheno_vals):
                y_values[col].append(val)
        except (ValueError, KeyError):
            continue

    n_samples = len(sample_ids)
    log_mem(f"Loaded {n_samples} samples with valid phenotypes")

    if covar_file:
        X = np.column_stack([np.ones(n_samples), np.array(X_list)])
    else:
        X = np.ones((n_samples, 1))

    return sample_ids, sample_aliases, y_values, X


def build_rif_residuals(y_values: dict[str, list[float] | np.ndarray], pheno_cols: list[str],
                        X: np.ndarray, taus: np.ndarray) -> tuple[list[np.ndarray], np.ndarray, list[np.ndarray]]:
    """Compute residualized RIF matrices for one or more phenotypes."""
    log_mem(f"Residualizing phenotypes on {X.shape[1]} covariates...")
    Q, _ = np.linalg.qr(X)

    rif_resid_list = []
    q_tau_list = []
    for col in pheno_cols:
        log_mem(f"Computing RIF matrix for {col}...")
        t_rif = time.perf_counter()
        y = np.asarray(y_values[col], dtype=np.float64)
        rif, q_tau, _ = compute_rif(y, taus)
        log_mem(f"RIF build for {col}: {time.perf_counter() - t_rif:.2f}s")
        rif_resid_list.append(rif - Q @ (Q.T @ rif))
        q_tau_list.append(q_tau)

    return rif_resid_list, Q, q_tau_list


class VariantGenerator:
    """Generator that yields variants one at a time while keeping bfile handle alive."""
    def __init__(self, bfile, snp_set):
        self.bfile = bfile  # Keep handle alive!
        self.snp_set = snp_set
        self.found = 0
        self.target = len(snp_set)
        
    def __iter__(self):
        for var in self.bfile:
            if var.rsid in self.snp_set:
                self.found += 1
                yield var
                if self.found >= self.target:
                    break


def get_variants_bgen(bgen_path: str, sample_path: Optional[str], 
                      snp_list: Optional[List[str]]):
    """
    Get variants from BGEN file using the bgen package.
    
    Uses delay_parsing=True for efficient streaming and pre-filters SNPs
    using the BGEN index to avoid scanning variants not in the file.
    
    Returns
    -------
    variants : iterable of variant objects
    samples : list of sample IDs
    bfile_handle : BgenReader handle (keep alive during iteration)
    """
    from bgen import BgenReader
    
    log_mem("Using bgen package with delay_parsing...")
    bfile = BgenReader(bgen_path, sample_path=sample_path or '', delay_parsing=True)
    
    if snp_list:
        log_mem(f"Fetching {len(snp_list)} specific SNPs using bgen...")
        
        # Pre-filter: intersect user SNPs with BGEN contents using index
        log_mem("Reading BGEN index for variant list...")
        bgen_rsids = set(bfile.rsids())
        snp_set = set(snp_list)
        valid_snps = snp_set & bgen_rsids
        invalid_count = len(snp_set) - len(valid_snps)
        
        if invalid_count > 0:
            log_mem(f"{invalid_count} requested SNPs not in BGEN, processing {len(valid_snps)}")
        
        if len(valid_snps) == 0:
            logger.warning("No valid SNPs found in BGEN file!")
        
        # Use generator for memory-efficient streaming
        gen = VariantGenerator(bfile, valid_snps)
        return gen, bfile.samples, bfile
    else:
        # Return the whole file for iteration
        log_mem("No SNP list provided; streaming full BGEN directly.")
        return bfile, bfile.samples, bfile


def match_samples(bgen_samples: list, sample_ids: list):
    """Match phenotype samples to BGEN samples and return indices."""
    sample_to_idx = {s: i for i, s in enumerate(bgen_samples)}
    matched_indices = [sample_to_idx.get(fid) for fid in sample_ids]
    valid_mask = np.array([idx is not None for idx in matched_indices])
    bgen_indices = np.array([idx for idx in matched_indices if idx is not None])
    return valid_mask, bgen_indices


def extract_dosage(variant, bgen_indices: Optional[np.ndarray] = None):
    """
    Extract dosage from variant for the first allele (alleles[0]).
    This ensures dosage matches effect_allele which is set to alleles[0].
    
    BGEN dosage is the expected count of allele[1] (alt_dosage).
    Dosage for allele[0] = 2 - alt_dosage (assuming diploid).
    
    Parameters
    ----------
    variant : variant object from bgen
    bgen_indices : indices to extract, or None for all samples
    
    Returns
    -------
    dosages : numpy array of dosages for alleles[0]
    """
    # alt_dosage is the expected count of allele[1]
    # For diploid, dosage of allele[0] = 2 - alt_dosage
    alt_dosages = variant.alt_dosage
    dosages = 2.0 - alt_dosages
    if bgen_indices is not None:
        return dosages[bgen_indices]
    return dosages


def write_stats_batch(stats_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
                      f_cov_scale_list, all_betas_list, n_cov_elements, n_taus):
    """Write one processed batch to the stage1 outputs."""
    for p_idx, stats in enumerate(stats_list):
        f_tsv = f_tsv_list[p_idx]
        f_cov = f_cov_list[p_idx]
        f_cov_ids = f_cov_ids_list[p_idx]
        f_cov_scale = f_cov_scale_list[p_idx]
        all_betas = all_betas_list[p_idx]

        for snp_stats, info in zip(stats, batch_info):
            if np.isnan(snp_stats[0, 0]):
                continue

            beta, se, Sigma = compute_covariance_from_blocks(snp_stats, args.blocks)

            row = list(info)
            for b, s in zip(beta, se):
                row.extend([f'{b:.6g}', f'{s:.6g}'])
            f_tsv.write('\t'.join(map(str, row)) + '\n')

            cov_upper = Sigma[np.triu_indices(n_taus)]
            if args.cov_dtype == 'int8':
                scale = np.max(np.abs(cov_upper)) / 127.0
                if not np.isfinite(scale) or scale <= 0:
                    scale = 1.0
                cov_q = np.clip(np.rint(cov_upper / scale), -127, 127).astype(np.int8)
                f_cov.write(cov_q.tobytes(order='C'))
                f_cov_scale.write(struct.pack('f', float(scale)))
            else:
                f_cov.write(struct.pack(f'{n_cov_elements}f', *cov_upper))
            f_cov_ids.write(f"{info[0]}\n")

            if args.output_rtau:
                all_betas.append(beta)


def compute_chr_rif_residuals(pheno_cols: list[str], y_values_full: dict[str, np.ndarray],
                              X_full: np.ndarray, taus: np.ndarray,
                              loco_predictions: Optional[dict[str, dict[str, np.ndarray]]],
                              chrom: Optional[str], valid_mask: np.ndarray) -> tuple[list[np.ndarray], np.ndarray, list[np.ndarray]]:
    """Build residualized RIF matrices for a chromosome, optionally using LOCO-adjusted phenotypes."""
    if loco_predictions is None:
        y_values_chr = {col: y_values_full[col][valid_mask] for col in pheno_cols}
    else:
        chrom_key = normalize_chr_label(chrom)
        y_values_chr = {}
        for col in pheno_cols:
            per_chr = loco_predictions[col]
            if chrom_key not in per_chr:
                raise ValueError(f"Chromosome {chrom_key} missing from LOCO predictions for phenotype '{col}'")
            y_values_chr[col] = y_values_full[col] - per_chr[chrom_key]
        y_values_chr = {col: y_values_chr[col][valid_mask] for col in pheno_cols}

    rif_resid_list, q_chr, q_tau_list = build_rif_residuals(y_values_chr, pheno_cols, X_full[valid_mask], taus)
    return rif_resid_list, q_chr, q_tau_list


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
                        help='Comma-separated phenotype columns to analyze in one pass '
                             '(uses only samples with non-missing values across all phenotypes/covariates)')
    parser.add_argument('--out', required=True,
                        help='Output prefix (creates {out}.stage1.tsv.gz, {out}.cov.gz)')
    
    # Optional arguments
    parser.add_argument('--sample', default=None,
                        help='Path to sample file (if not embedded in BGEN)')
    parser.add_argument('--covar', default=None,
                        help='Path to covariate file')
    parser.add_argument('--covar-col', dest='covar_cols', nargs='+', default=None,
                        help='Covariate column(s) to use (space or comma-separated). '
                             'If omitted, all non-ID columns in the covariate file are used.')
    parser.add_argument('--snps', default=None,
                        help='File with SNP IDs to process (one per line)')
    parser.add_argument('--keep', default=None,
                        help='File with sample IDs to keep (one per line)')
    parser.add_argument('--loco-preds', dest='loco_preds', default=None,
                        help='REGENIE step1 LOCO predictions (for future use)')
    parser.add_argument('--allow-numpy-fallback', action='store_true',
                        help='Allow the uncompiled NumPy fallback. This can be ~1000x slower and is for debugging only.')
    
    # Algorithm parameters
    default_taus = ",".join([f"{t:.3f}" for t in np.linspace(0.02, 0.98, 45)])
    parser.add_argument('--taus', default=default_taus,
                        help='Comma-separated tau values (default: 45 taus from 0.02 to 0.98)')
    parser.add_argument('--blocks', type=int, default=None,
                        help='Number of jackknife blocks (default: max(taus + 20, 64))')
    parser.add_argument('--batch-size', dest='batch_size', type=int, default=500,
                        help='SNPs per batch for memory efficiency')
    parser.add_argument('--threads', type=int, default=1,
                        help='OpenMP threads for C++ kernel')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for block assignment')
    
    # Output options
    parser.add_argument('--output-rtau', dest='output_rtau', action='store_true',
                        help='Also output R_tau correlation matrix for plugin_cor')
    parser.add_argument('--cov-dtype', dest='cov_dtype', choices=['int8', 'float32'],
                        default='int8',
                        help='Covariance storage dtype: int8 (with .cov.scale.gz sidecar) '
                             'or float32 (legacy). Default: int8')
    
    args = parser.parse_args()
    set_numpy_fallback_allowed(args.allow_numpy_fallback)
    require_cpp_extension("fungwas-stage1 CLI")

    if not args.pheno_col and not args.pheno_cols:
        parser.error("Provide --pheno-col or --pheno-cols")
    
    # Parse taus
    taus = np.array([float(t) for t in args.taus.split(',')])
    T = len(taus)
    if args.blocks is None:
        args.blocks = max(T + 20, 64)
        logger.info("Using auto --blocks=%d (max(taus + 20, 64))", args.blocks)

    if args.blocks <= T:
        parser.error(
            f"--blocks ({args.blocks}) must be greater than number of taus ({T}). "
            "Increase --blocks or reduce --taus."
        )
    if args.blocks < T + 20:
        logger.warning(
            "--blocks=%d is only slightly above taus=%d; covariance may be noisy. "
            "Consider --blocks >= taus + 20.",
            args.blocks, T
        )
    
    log_mem(f"Starting FungWas Stage 1 (taus={T}, blocks={args.blocks})")
    logger.info("Note: For Stage 2, consider mean-centering non-intercept columns of your W matrix to reduce collinearity.")
    t_start = time.time()
    t_load = time.perf_counter()
    
    if args.pheno_cols:
        pheno_cols = [c.strip() for c in args.pheno_cols.split(',') if c.strip()]
    else:
        pheno_cols = [args.pheno_col]
    covar_cols = parse_list_arg(args.covar_cols)

    use_multi = len(pheno_cols) > 1
    if not use_multi and args.pheno_cols:
        log_mem("Single phenotype provided via --pheno-cols; using single-phenotype path.")

    if use_multi:
        log_mem(f"Loading {len(pheno_cols)} phenotypes in one pass...")
    sample_ids, sample_aliases, y_values_raw, X_full = load_analysis_inputs(
        args.pheno, pheno_cols, args.covar, keep_file=args.keep, covar_cols=covar_cols
    )
    y_values_full = {col: np.asarray(values, dtype=np.float64) for col, values in y_values_raw.items()}
    loco_predictions = None
    if args.loco_preds:
        loco_predictions, loco_valid_mask = load_loco_predictions(args.loco_preds, pheno_cols, sample_aliases)
        if not np.all(loco_valid_mask):
            dropped = int((~loco_valid_mask).sum())
            log_mem(f"Dropping {dropped} phenotype rows not present in LOCO predictions")
            sample_ids = [sid for sid, keep in zip(sample_ids, loco_valid_mask) if keep]
            sample_aliases = [aliases for aliases, keep in zip(sample_aliases, loco_valid_mask) if keep]
            y_values_full = {col: values[loco_valid_mask] for col, values in y_values_full.items()}
            X_full = X_full[loco_valid_mask]
        log_mem("LOCO predictions loaded; RIFs will be recomputed per chromosome on Y - LOCO(chr)")
    RIF_resid_list, Q, q_tau_list = build_rif_residuals(y_values_full, pheno_cols, X_full, taus)
    log_mem(f"Phenotype load time: {time.perf_counter() - t_load:.2f}s")
    N = len(sample_ids)
    
    # Setup output dir
    out_dir = os.path.dirname(args.out) or '.'
    os.makedirs(out_dir, exist_ok=True)
    
    # Load SNP list if provided
    snp_list = None
    if args.snps:
        snp_list = load_snps_list(args.snps)
        log_mem(f"Loaded {len(snp_list)} SNPs from {args.snps}")
    
    # Open BGEN and get variants
    variants, bgen_samples, bfile_handle = get_variants_bgen(
        args.bgen, args.sample, snp_list
    )
    
    # Match samples
    valid_mask, bgen_indices = match_samples(bgen_samples, sample_ids)
    
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
    out_cov_scale_list = [f"{prefix}.cov.scale.gz" for prefix in out_prefixes]
    out_cov_meta_list = [f"{prefix}.cov.meta.json" for prefix in out_prefixes]
    out_cov_ids_list = [f"{prefix}.cov.ids.tsv.gz" for prefix in out_prefixes]
    
    # Column headers
    cols = ['snp_id', 'chr', 'bp', 'effect_allele', 'other_allele']
    for t in taus:
        # Use 3 decimals to avoid tau label collisions on dense grids (e.g. 100 taus).
        cols.extend([f'beta_tau_{t:.3f}', f'se_tau_{t:.3f}'])
    
    # Upper triangle size
    n_cov_elements = T * (T + 1) // 2
    
    # For R_tau computation
    all_betas_list = [([] if args.output_rtau else None) for _ in pheno_cols]
    
    # Process in batches
    total_snps = 0
    skipped_non_biallelic = 0
    G_batch = np.zeros((N_matched, args.batch_size), dtype=np.float64, order='F')
    batch_info = []
    
    from contextlib import ExitStack

    total_compute = 0.0
    total_write = 0.0
    active_chr = None

    with ExitStack() as stack:
        f_tsv_list = [stack.enter_context(gzip.open(path, 'wt')) for path in out_tsv_list]
        f_cov_list = [stack.enter_context(gzip.open(path, 'wb')) for path in out_cov_list]
        f_cov_scale_list = (
            [stack.enter_context(gzip.open(path, 'wb')) for path in out_cov_scale_list]
            if args.cov_dtype == 'int8' else [None for _ in out_cov_scale_list]
        )
        f_cov_ids_list = [stack.enter_context(gzip.open(path, 'wt')) for path in out_cov_ids_list]
        for f_tsv in f_tsv_list:
            f_tsv.write('\t'.join(cols) + '\n')
        for f_ids in f_cov_ids_list:
            f_ids.write('snp_id\n')
        
        idx_in_batch = 0
        
        for variant in variants:
            if len(variant.alleles) != 2:
                skipped_non_biallelic += 1
                continue

            variant_chr = normalize_chr_label(variant.chrom)
            if loco_predictions is not None and variant_chr != active_chr:
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
                    write_stats_batch(
                        stats_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
                        f_cov_scale_list, all_betas_list, n_cov_elements, T
                    )
                    total_write += time.perf_counter() - t_write
                    total_snps += idx_in_batch
                    idx_in_batch = 0
                    batch_info = []

                log_mem(f"Applying LOCO phenotype adjustment for chromosome {variant_chr}")
                RIF_resid_list, Q, q_tau_list = compute_chr_rif_residuals(
                    pheno_cols, y_values_full, X_full, taus, loco_predictions, variant_chr, valid_mask
                )
                active_chr = variant_chr

            # Extract dosage
            all_dosages = extract_dosage(variant)
            dosages = all_dosages[bgen_indices]
            
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
                write_stats_batch(
                    stats_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
                    f_cov_scale_list, all_betas_list, n_cov_elements, T
                )
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
            write_stats_batch(
                stats_list, batch_info, args, f_tsv_list, f_cov_list, f_cov_ids_list,
                f_cov_scale_list, all_betas_list, n_cov_elements, T
            )
            total_write += time.perf_counter() - t_write
            
            total_snps += idx_in_batch
    
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
                f.write('\t'.join([f'tau_{t:.3f}' for t in taus]) + '\n')
                for row in R_tau:
                    f.write('\t'.join([f'{v:.6f}' for v in row]) + '\n')
            logger.info(f"Wrote R_tau: {rtau_file}")
    
    elapsed = time.time() - t_start
    for meta_path in out_cov_meta_list:
        with open(meta_path, 'w') as f:
            json.dump({
                "cov_dtype": args.cov_dtype,
                "n_taus": int(T),
                "n_cov_elements": int(n_cov_elements),
                "record_bytes": int(n_cov_elements if args.cov_dtype == 'int8' else n_cov_elements * 4),
                "scale_file": bool(args.cov_dtype == 'int8')
            }, f, indent=2)
    log_mem(
        f"Batch compute time: {total_compute:.1f}s, write time: {total_write:.1f}s, "
        f"batch total: {time.perf_counter() - t_batches:.1f}s"
    )
    log_mem(f"Dropped {skipped_non_biallelic} non-biallelic variants")
    log_mem(f"Done! {total_snps} SNPs in {elapsed:.1f}s ({total_snps/elapsed:.1f} SNPs/s)")
    logger.info(f"Output: {', '.join(out_tsv_list)}")


if __name__ == "__main__":
    main()
