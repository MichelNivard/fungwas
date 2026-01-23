"""
Simulate data and run stage1 (Python) to generate Stage 1 outputs.

Creates:
- {out}.stage1.tsv.gz with beta/SE per tau
- {out}.cov.gz with per-SNP covariance upper triangles (float32)
- {out}.meta.json with simulation parameters and true effects
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from fungwas_stage1.core import run_stage1


def simulate_data(
    n_samples: int,
    n_snps: int,
    n_covars: int,
    seed: int,
):
    rng = np.random.default_rng(seed)

    maf = rng.uniform(0.05, 0.5, size=n_snps)
    G = rng.binomial(2, maf, size=(n_samples, n_snps)).astype(float)

    covariates = rng.normal(size=(n_samples, n_covars)) if n_covars > 0 else None

    p1 = 0.6
    mu1, sd1 = 0.0, 0.8
    mu2, sd2 = 3.0, 1.2

    beta_gamma_true = rng.normal(size=n_snps) / 30
    beta_1_true = rng.normal(size=n_snps) / 30
    beta_2_true = rng.normal(size=n_snps) / 30

    cov_effect = 0.0
    cov_beta = None
    if covariates is not None:
        cov_beta = rng.normal(size=n_covars) / 4
        cov_effect = covariates @ cov_beta

    gamma = np.log(p1 / (1 - p1)) + G @ beta_gamma_true
    pi_i = 1 / (1 + np.exp(-gamma))
    z = rng.binomial(1, pi_i, size=n_samples)

    Y = np.empty(n_samples)
    idx1 = z == 1
    idx2 = ~idx1

    Y[idx1] = rng.normal(
        loc=mu1 + G[idx1] @ beta_1_true + cov_effect[idx1],
        scale=sd1,
    )
    Y[idx2] = rng.normal(
        loc=mu2 + G[idx2] @ beta_2_true + cov_effect[idx2],
        scale=sd2,
    )

    params = {
        "p1": p1,
        "mu1": mu1,
        "sd1": sd1,
        "mu2": mu2,
        "sd2": sd2,
        "beta_gamma_true": beta_gamma_true.tolist(),
        "beta_1_true": beta_1_true.tolist(),
        "beta_2_true": beta_2_true.tolist(),
        "cov_beta": None if cov_beta is None else cov_beta.tolist(),
    }

    return G, Y, covariates, params


def write_stage1_outputs(
    out_prefix: Path,
    taus: np.ndarray,
    betas: np.ndarray,
    ses: np.ndarray,
    covariances: list[np.ndarray],
):
    out_tsv = out_prefix.with_suffix(".stage1.tsv.gz")
    out_cov = out_prefix.with_suffix(".cov.gz")

    columns = ["snp_id", "chr", "bp", "effect_allele", "other_allele"]
    for tau in taus:
        columns.append(f"beta_tau_{tau:.2f}")
        columns.append(f"se_tau_{tau:.2f}")

    n_snps, n_taus = betas.shape
    triu_idx = np.triu_indices(n_taus)

    with gzip.open(out_tsv, "wt") as f_tsv, gzip.open(out_cov, "wb") as f_cov:
        f_tsv.write("\t".join(columns) + "\n")

        for j in range(n_snps):
            row = [f"SNP{j+1}", 1, j + 1, "A", "G"]
            for b, s in zip(betas[j], ses[j]):
                row.extend([f"{b:.6g}", f"{s:.6g}"])
            f_tsv.write("\t".join(map(str, row)) + "\n")

            cov_upper = covariances[j][triu_idx].astype(np.float32)
            f_cov.write(cov_upper.tobytes())

    return out_tsv, out_cov


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--n-samples", type=int, default=300000)
    parser.add_argument("--n-snps", type=int, default=100)
    parser.add_argument("--n-covars", type=int, default=2)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--taus", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9")
    parser.add_argument("--n-blocks", type=int, default=25)
    parser.add_argument("--n-threads", type=int, default=4)
    args = parser.parse_args()

    taus = np.array([float(x) for x in args.taus.split(",")])

    G, Y, covariates, params = simulate_data(
        args.n_samples, args.n_snps, args.n_covars, args.seed
    )

    result = run_stage1(
        G,
        Y,
        taus,
        covariates=covariates,
        n_blocks=args.n_blocks,
        seed=args.seed,
        n_threads=args.n_threads,
    )

    out_prefix = Path(args.out_prefix)
    out_tsv, out_cov = write_stage1_outputs(
        out_prefix,
        taus,
        result["betas"],
        result["se"],
        result["covariances"],
    )

    meta = {
        "taus": taus.tolist(),
        "q_tau": result["q_tau"].tolist(),
        **params,
    }

    meta_path = out_prefix.with_suffix(".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2))

    print(f"Wrote stage1: {out_tsv}")
    print(f"Wrote cov: {out_cov}")
    print(f"Wrote meta: {meta_path}")


if __name__ == "__main__":
    main()
