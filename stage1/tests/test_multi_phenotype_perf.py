import time

import numpy as np

from fungwas_stage1 import _stage1_cpp


def rss_mb() -> float:
    with open("/proc/self/statm") as handle:
        parts = handle.read().strip().split()
    rss_pages = int(parts[1])
    page_size = 4096
    return rss_pages * page_size / (1024 * 1024)


def run_case(num_pheno: int, n_threads: int = 2):
    print(f"\n--- Running case: {num_pheno} phenotype(s) ---", flush=True)
    total_start = time.perf_counter()
    rng = np.random.default_rng(123)
    n, m, t = 2000, 2000, 5
    maf = rng.uniform(0.1, 0.4, size=m)
    g = rng.binomial(2, maf, size=(n, m)).astype(np.float64)
    g[rng.random(size=g.shape) < 0.02] = np.nan
    print(f"Data generation: {time.perf_counter() - total_start:.3f}s", flush=True)

    x = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
    q, _ = np.linalg.qr(x)
    print(f"Covariate prep: {time.perf_counter() - total_start:.3f}s", flush=True)

    target_snp = 0
    true_beta = 0.2
    g_target = g[:, target_snp].copy()
    mean_target = np.nanmean(g_target)
    g_target[np.isnan(g_target)] = mean_target
    rif_list = []
    for _ in range(num_pheno):
        noise = rng.normal(scale=0.5, size=(n, t))
        rif = g_target[:, None] * true_beta + noise
        rif_list.append(np.asfortranarray(rif, dtype=np.float64))
    print(f"RIF prep: {time.perf_counter() - total_start:.3f}s", flush=True)

    g_f = np.asfortranarray(g)
    q_f = np.asfortranarray(q)

    rss_before = rss_mb()
    start = time.perf_counter()
    results = _stage1_cpp.process_block_dense_impl(
        g_f,
        rif_list,
        q_f,
        5.0,
        1e-8,
        n_threads,
    )
    elapsed = time.perf_counter() - start
    rss_after = rss_mb()
    print(f"C++ kernel time: {elapsed:.3f}s", flush=True)
    print(f"RSS delta: {rss_after - rss_before:.3f} MB", flush=True)

    betas = results[0][0]
    approx = np.nanmean(np.abs(betas[target_snp, :] - true_beta))

    return {
        "elapsed_sec": elapsed,
        "rss_before_mb": rss_before,
        "rss_after_mb": rss_after,
        "rss_delta_mb": rss_after - rss_before,
        "mean_abs_error": approx,
        "total_sec": time.perf_counter() - total_start,
    }


if __name__ == "__main__":
    one = run_case(1)
    two = run_case(2)

    print("1 phenotype:", one)
    print("2 phenotypes:", two)
    print("extra_rss_mb:", two["rss_after_mb"] - one["rss_after_mb"])
