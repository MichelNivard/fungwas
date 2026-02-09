#!/usr/bin/env python3
"""
Simulate full raw-data Stage1 input with known per-SNP model truth.

This script simulates:
- Genotypes G (optionally with LD blocks)
- Phenotype Y from a quantile-effect model where each SNP has a known truth model
  among the supplied weight matrices (e.g. add/shift/both/vqtl)

Then it runs real Stage1 (`fungwas_stage1.core.run_stage1`) and writes:
- {out_prefix}.stage1.tsv.gz
- {out_prefix}.cov.gz
- {out_prefix}.cov.ids.tsv.gz
- {out_prefix}.truth.tsv.gz         (true model, true theta, true tau-slope profile)
- {out_prefix}.meta.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from fungwas_stage1.core import run_stage1


def parse_weights_arg(s: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"Malformed weights entry: {part}")
        name, path = part.split("=", 1)
        name = name.strip()
        path = path.strip()
        if not name or not path:
            raise ValueError(f"Malformed weights entry: {part}")
        out[name] = path
    if not out:
        raise ValueError("No weights parsed from --weights")
    return out


def parse_model_probs(s: str, model_names: List[str]) -> Dict[str, float]:
    probs = {m: 0.0 for m in model_names}
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"Malformed model prob entry: {part}")
        name, val = part.split(":", 1)
        name = name.strip()
        if name not in probs:
            raise ValueError(f"Unknown model in --model-probs: {name}")
        probs[name] = float(val)
    total = sum(probs.values())
    if total <= 0:
        raise ValueError("Model probabilities sum to zero")
    return {k: v / total for k, v in probs.items()}


def load_weight(path: str) -> Tuple[np.ndarray, np.ndarray]:
    # RDS files are read via pyreadr when available.
    # Fallback: call Rscript to emit a plain-text matrix with tau first column.
    try:
        import pyreadr
        res = pyreadr.read_r(path)
        if not res:
            raise ValueError(f"Unable to read weight file: {path}")
        obj = next(iter(res.values()))

        W = None
        taus = None
        if isinstance(obj, dict):
            if "W" in obj:
                W = np.asarray(obj["W"], dtype=float)
            if "taus" in obj:
                taus = np.asarray(obj["taus"], dtype=float)
        else:
            W = np.asarray(obj, dtype=float)

        if W is None:
            raise ValueError(f"Could not parse W from {path}")
        if taus is not None and W.shape[0] != len(taus) and W.shape[1] == len(taus):
            W = W.T
        if taus is None:
            raise ValueError(f"Weight file missing taus: {path}")
        return W, taus
    except ImportError:
        with tempfile.NamedTemporaryFile(prefix="weight_", suffix=".tsv", delete=False) as tf:
            tmp_path = tf.name
        r_code = (
            "x <- readRDS(commandArgs(TRUE)[1]); "
            "out <- commandArgs(TRUE)[2]; "
            "W <- if (is.list(x) && !is.null(x$W)) x$W else x; "
            "taus <- if (is.list(x) && !is.null(x$taus)) x$taus else NULL; "
            "if (is.null(taus)) stop('weights must include taus'); "
            "if (nrow(W) != length(taus) && ncol(W) == length(taus)) W <- t(W); "
            "M <- cbind(tau=as.numeric(taus), as.matrix(W)); "
            "write.table(M, file=out, sep='\\t', row.names=FALSE, col.names=FALSE, quote=FALSE)"
        )
        try:
            subprocess.run(
                ["Rscript", "-e", r_code, path, tmp_path],
                check=True,
                capture_output=True,
                text=True,
            )
            mat = np.loadtxt(tmp_path, delimiter="\t")
            if mat.ndim == 1:
                mat = mat.reshape(1, -1)
            taus = np.asarray(mat[:, 0], dtype=float)
            W = np.asarray(mat[:, 1:], dtype=float)
            return W, taus
        finally:
            try:
                Path(tmp_path).unlink(missing_ok=True)
            except Exception:
                pass


def align_weights(weight_paths: Dict[str, str]) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    loaded = {k: load_weight(v) for k, v in weight_paths.items()}
    first = next(iter(loaded))
    taus_ref = np.asarray(loaded[first][1], dtype=float)
    ord_idx = np.argsort(taus_ref)
    taus_ref = taus_ref[ord_idx]

    aligned: Dict[str, np.ndarray] = {}
    for name, (W, taus) in loaded.items():
        taus = np.asarray(taus, dtype=float)
        if W.shape[0] != len(taus) and W.shape[1] == len(taus):
            W = W.T

        idx = []
        for t in taus_ref:
            j = int(np.argmin(np.abs(taus - t)))
            if abs(taus[j] - t) > 1e-8:
                raise ValueError(f"Tau mismatch for model {name} at tau={t}")
            idx.append(j)
        aligned[name] = np.asarray(W[idx, :], dtype=float)

    return taus_ref, aligned


def gen_ld_genotypes(
    n_samples: int,
    n_snps: int,
    rng: np.random.Generator,
    ld_block: int,
    ld_rho: float,
    maf_min: float,
    maf_max: float,
) -> np.ndarray:
    maf = rng.uniform(maf_min, maf_max, size=n_snps)
    G = np.zeros((n_samples, n_snps), dtype=float)

    for start in range(0, n_snps, ld_block):
        end = min(start + ld_block, n_snps)
        b = end - start
        if b == 1 or ld_rho <= 0:
            z = rng.normal(size=(n_samples, b))
        else:
            idx = np.arange(b)
            R = ld_rho ** np.abs(np.subtract.outer(idx, idx))
            L = np.linalg.cholesky(R)
            z = rng.normal(size=(n_samples, b)) @ L.T

        for j in range(b):
            p = maf[start + j]
            p0 = (1 - p) ** 2
            p1 = 2 * p * (1 - p)
            t0 = np.quantile(z[:, j], p0)
            t1 = np.quantile(z[:, j], p0 + p1)
            g = np.where(z[:, j] < t0, 0.0, np.where(z[:, j] < t1, 1.0, 2.0))
            G[:, start + j] = g

    return G


def linear_interp_basis(u: np.ndarray, taus: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Returns left index, right index, right weight.
    u = np.clip(u, taus[0], taus[-1])
    r = np.searchsorted(taus, u, side="left")
    r = np.clip(r, 1, len(taus) - 1)
    l = r - 1
    denom = taus[r] - taus[l]
    wr = (u - taus[l]) / denom
    return l, r, wr


def write_stage1_outputs(
    out_prefix: Path,
    taus: np.ndarray,
    betas: np.ndarray,
    ses: np.ndarray,
    covariances: List[np.ndarray],
) -> Tuple[Path, Path, Path]:
    out_tsv = out_prefix.with_suffix(".stage1.tsv.gz")
    out_cov = out_prefix.with_suffix(".cov.gz")
    out_cov_ids = out_prefix.with_suffix(".cov.ids.tsv.gz")

    columns = ["snp_id", "chr", "bp", "effect_allele", "other_allele"]
    for tau in taus:
        columns.append(f"beta_tau_{tau:.6g}")
        columns.append(f"se_tau_{tau:.6g}")

    n_snps, n_taus = betas.shape
    triu_idx = np.triu_indices(n_taus)

    with gzip.open(out_tsv, "wt") as f_tsv, gzip.open(out_cov, "wb") as f_cov, gzip.open(out_cov_ids, "wt") as f_ids:
        f_tsv.write("\t".join(columns) + "\n")
        f_ids.write("snp_id\tchr\tbp\teffect_allele\tother_allele\n")

        for j in range(n_snps):
            sid = f"SNP{j+1}"
            row = [sid, 1, j + 1, "A", "G"]
            for b, s in zip(betas[j], ses[j]):
                row.extend([f"{b:.9g}", f"{s:.9g}"])
            f_tsv.write("\t".join(map(str, row)) + "\n")
            f_ids.write(f"{sid}\t1\t{j+1}\tA\tG\n")

            cov_upper = np.asarray(covariances[j])[triu_idx].astype(np.float32)
            f_cov.write(cov_upper.tobytes())

    return out_tsv, out_cov, out_cov_ids


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--weights", required=True, help="name=path.rds,name=path.rds,...")
    ap.add_argument("--model-probs", default="add:0.25,shift:0.25,both:0.25,vqtl:0.25")
    ap.add_argument("--n-samples", type=int, default=12000)
    ap.add_argument("--n-snps", type=int, default=1200)
    ap.add_argument("--n-covars", type=int, default=2)
    ap.add_argument("--ld-block", type=int, default=20)
    ap.add_argument("--ld-rho", type=float, default=0.2)
    ap.add_argument("--maf-min", type=float, default=0.05)
    ap.add_argument("--maf-max", type=float, default=0.5)
    ap.add_argument("--theta-sd", type=float, default=0.05)
    ap.add_argument("--baseline", choices=["normal", "lognormal"], default="normal")
    ap.add_argument("--resid-sd", type=float, default=0.05)
    ap.add_argument("--n-taus", type=int, default=45,
                    help="Number of taus to use from the aligned weight grid (default: 45)")
    ap.add_argument("--tau-min", type=float, default=0.02,
                    help="Minimum tau to include before subsetting (default: 0.02)")
    ap.add_argument("--tau-max", type=float, default=0.98,
                    help="Maximum tau to include before subsetting (default: 0.98)")
    ap.add_argument("--n-blocks", type=int, default=-1,
                    help="Jackknife blocks; default auto max(n_taus + 20, 64)")
    ap.add_argument("--n-threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=6520)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    weight_paths = parse_weights_arg(args.weights)
    taus, W_list = align_weights(weight_paths)
    elig = np.where((taus >= args.tau_min) & (taus <= args.tau_max))[0]
    if len(elig) < 2:
        raise ValueError("Tau filter leaves fewer than 2 taus; adjust --tau-min/--tau-max")
    n_taus = int(min(max(2, args.n_taus), len(elig)))
    pick_loc = np.unique(np.round(np.linspace(0, len(elig) - 1, n_taus)).astype(int))
    pick = elig[pick_loc]
    taus = taus[pick]
    W_list = {k: v[pick, :] for k, v in W_list.items()}
    model_names = list(W_list.keys())
    model_probs = parse_model_probs(args.model_probs, model_names)

    n_samples = args.n_samples
    n_snps = args.n_snps
    if args.n_blocks <= 0:
        args.n_blocks = max(len(taus) + 20, 64)
    if args.n_blocks <= len(taus):
        raise ValueError(
            f"n_blocks ({args.n_blocks}) must be > number of taus ({len(taus)}). "
            "Use more blocks or fewer taus."
        )

    G = gen_ld_genotypes(
        n_samples=n_samples,
        n_snps=n_snps,
        rng=rng,
        ld_block=max(1, args.ld_block),
        ld_rho=max(0.0, min(0.99, args.ld_rho)),
        maf_min=args.maf_min,
        maf_max=args.maf_max,
    )

    covariates = rng.normal(size=(n_samples, args.n_covars)) if args.n_covars > 0 else None

    # Assign true model and theta per SNP.
    probs = np.array([model_probs[m] for m in model_names], dtype=float)
    probs = probs / probs.sum()
    model_idx = rng.choice(len(model_names), size=n_snps, replace=True, p=probs)
    true_model = np.array([model_names[i] for i in model_idx], dtype=object)

    max_k = max(W.shape[1] for W in W_list.values())
    theta_true = np.full((n_snps, max_k), np.nan, dtype=float)
    s_true = np.zeros((n_snps, len(taus)), dtype=float)

    for j in range(n_snps):
        m = true_model[j]
        W = W_list[m]
        k = W.shape[1]
        th = rng.normal(0.0, args.theta_sd, size=k)
        theta_true[j, :k] = th
        s_true[j, :] = W @ th

    # Build phenotype via baseline quantile + sum_j g_ij * s_j(u_i).
    u = rng.uniform(taus[0], taus[-1], size=n_samples)
    if args.baseline == "normal":
        y_base = rng.normal(0.0, 1.0, size=n_samples)
    else:
        y_base = np.exp(rng.normal(0.0, 0.5, size=n_samples))

    l_idx, r_idx, w_r = linear_interp_basis(u, taus)
    y = y_base.copy()

    # Add SNP contributions.
    for j in range(n_snps):
        s = s_true[j]
        f_u = s[l_idx] * (1.0 - w_r) + s[r_idx] * w_r
        y += G[:, j] * f_u

    if args.resid_sd > 0:
        y += rng.normal(0.0, args.resid_sd, size=n_samples)

    result = run_stage1(
        G,
        y,
        taus,
        covariates=covariates,
        n_blocks=args.n_blocks,
        seed=args.seed,
        n_threads=args.n_threads,
    )

    out_prefix = Path(args.out_prefix)
    out_tsv, out_cov, out_cov_ids = write_stage1_outputs(
        out_prefix,
        taus,
        result["betas"],
        result["se"],
        result["covariances"],
    )

    truth_path = out_prefix.with_suffix(".truth.tsv.gz")
    with gzip.open(truth_path, "wt") as f:
        cols = ["snp_id", "true_model"]
        cols.extend([f"theta_{k+1}" for k in range(max_k)])
        cols.extend([f"true_beta_tau_{t:.6g}" for t in taus])
        f.write("\t".join(cols) + "\n")

        for j in range(n_snps):
            row = [f"SNP{j+1}", true_model[j]]
            row.extend(["" if np.isnan(v) else f"{v:.9g}" for v in theta_true[j]])
            row.extend([f"{v:.9g}" for v in s_true[j]])
            f.write("\t".join(row) + "\n")

    meta = {
        "seed": args.seed,
        "n_samples": n_samples,
        "n_snps": n_snps,
        "n_covars": args.n_covars,
        "ld_block": args.ld_block,
        "ld_rho": args.ld_rho,
        "maf_min": args.maf_min,
        "maf_max": args.maf_max,
        "theta_sd": args.theta_sd,
        "baseline": args.baseline,
        "resid_sd": args.resid_sd,
        "n_blocks": args.n_blocks,
        "n_threads": args.n_threads,
        "model_probs": model_probs,
        "taus": taus.tolist(),
        "q_tau_stage1": result["q_tau"].tolist(),
        "weights": weight_paths,
        "n_taus_used": int(len(taus)),
        "tau_min": float(args.tau_min),
        "tau_max": float(args.tau_max),
    }
    meta_path = out_prefix.with_suffix(".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2))

    print(f"Wrote stage1: {out_tsv}")
    print(f"Wrote cov: {out_cov}")
    print(f"Wrote cov ids: {out_cov_ids}")
    print(f"Wrote truth: {truth_path}")
    print(f"Wrote meta: {meta_path}")


if __name__ == "__main__":
    main()
