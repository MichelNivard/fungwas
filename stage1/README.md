# FungWas Stage 1: Fast RIF Quantile GWAS

Python/C++ implementation for blazing-fast RIF quantile regression with jackknife covariance estimation.

## Installation

### Prerequisites

Stage 1 requires a **C++ extension** for fast jackknife covariance computation. Without it, the numpy fallback is ~100x slower for large datasets.

```bash
# With conda (recommended)
conda env create -f environment.yml
conda activate fungwas-stage1

# Build C++ extension (REQUIRED for performance)
python setup.py build_ext --inplace

# Install package
pip install -e .
```

**Important:** If you see `C++ extension not found, using numpy fallback (slower)` when running, the C++ extension is not built. Build it with:
```bash
cd stage1
python setup.py build_ext --inplace
```

### Conda environment notes (BGEN support)

Stage 1 uses the Python `bgen` module for BGEN streaming, installed via pip (currently pinned to `bgen==1.9.7`). No external tools are required.

```bash
conda env update -n fungwas-stage1 -f environment.yml
```

## Usage

### Option 1: Array-based API (for testing / small data)

Similar to R's `quantile_gwas(Y, G, taus)`:

```python
import numpy as np
from fungwas_stage1 import run_stage1

# Simulate data
N, M = 5000, 100
G = np.random.binomial(2, 0.3, (N, M)).astype(float)
Y = np.random.randn(N)  # Raw phenotype (NOT pre-computed RIF!)
taus = np.linspace(0.02, 0.98, 45)

# Run Stage 1 - RIF is computed internally
result = run_stage1(G, Y, taus, n_blocks=None, n_threads=4)  # auto blocks = max(T+20, 64)

# Access results
print(result['betas'].shape)       # (M, T) tau-level betas
print(result['se'].shape)          # (M, T) standard errors
print(result['covariances'][0])    # (T, T) covariance for SNP 0
```

### Option 2: CLI for BGEN files (for HPC production)

```bash
# Basic usage - provide raw phenotype, RIF computed internally
fungwas-stage1 \
    --bgen data.chr22.bgen \
    --pheno phenotypes.txt \
    --pheno-col testosterone \
    --out results/chr22

# With all options
fungwas-stage1 \
    --bgen data.chr22.bgen \
    --sample data.sample \
    --pheno phenotypes.txt \
    --pheno-col testosterone \
    --covar covariates.txt \
    --snps chr22_hm3.txt \
    --taus "0.020,0.042,0.064,0.085,0.107,0.129,0.151,0.173,0.195,0.216,0.238,0.260,0.282,0.304,0.325,0.347,0.369,0.391,0.413,0.435,0.456,0.478,0.500,0.522,0.544,0.565,0.587,0.609,0.631,0.653,0.675,0.696,0.718,0.740,0.762,0.784,0.805,0.827,0.849,0.871,0.893,0.915,0.936,0.958,0.980" \
    --blocks 64 \
    --batch-size 500 \
    --threads 8 \
    --cov-dtype int8 \
    --output-rtau \
    --out results/chr22
```

#### Multi-phenotype CLI (single BGEN pass)

Provide comma-separated phenotype columns to scan multiple phenotypes with one
genotype pass. Outputs are suffixed with the phenotype name.

Note: the multi-phenotype path uses the intersection of samples with
non-missing values across all phenotypes (and covariates, if provided). This
means phenotypes with more missingness reduce the sample set for *all*
phenotypes in the joint run. This tradeoff enables shared residualization and
block computations for speed. If you need per-phenotype sample sizes, run each
phenotype separately using the single-phenotype CLI.

```bash
fungwas-stage1 \
    --bgen data.chr22.bgen \
    --pheno phenotypes.txt \
    --pheno-cols "phenotype_a,phenotype_b" \
    --covar covariates.txt \
    --snps chr22_hm3.txt \
    --taus "0.10,0.50,0.90" \
    --blocks 64 \
    --batch-size 500 \
    --threads 8 \
    --out results/chr22
```

If `--pheno-cols` contains a single phenotype, Stage 1 automatically falls back
to the single-phenotype kernel to avoid extra overhead.

### BGEN Loading Strategy

Stage 1 uses an efficient hybrid approach for BGEN loading:

1. **Pre-filtering**: Uses the BGEN index (`bfile.rsids()`) to intersect the user's SNP list with variants actually in the BGEN file
2. **delay_parsing**: Opens BGEN with `delay_parsing=True` to defer genotype decompression until needed
3. **Streaming**: Iterates through the BGEN file, yielding only matching variants
4. **Batching**: Processes variants in chunks (default 500) for memory efficiency

**Performance:**
- ~800 SNPs/second for full chromosome scans
- ~22 seconds for 17k HapMap3 SNPs on chr22 (463k samples)
- Constant memory usage (~140MB regardless of SNP count)

```bash
# Full chromosome (streams all variants)
fungwas-stage1 --bgen data.chr22.bgen --pheno pheno.txt --pheno-col y --out results/chr22

# SNP subset (pre-filters using index, then streams matches)
fungwas-stage1 --bgen data.chr22.bgen --pheno pheno.txt --pheno-col y --snps hm3_snps.txt --out results/chr22
```

### Multi-phenotype C++ kernel (advanced)

The C++ pybind module exposes `process_block_dense_impl`/`process_block_sparse_impl`
for running multiple RIF matrices against the same genotype block. This shares
the genotype residualization and leverage calculations across phenotypes.

```python
from fungwas_stage1 import _stage1_cpp

results = _stage1_cpp.process_block_dense_impl(G, [rif1, rif2], Q)
for betas, ses in results:
    print(betas.shape, ses.shape)
```

## Input Files

### Phenotype file (space/tab delimited)

Standard format with FID, IID, and raw phenotype column(s):
```
FID IID testosterone age_at_assessment
1001 1001 15.3 45
1002 1002 12.1 52
```

The `--pheno-col` argument specifies which column to analyze. FungWas computes RIF internally.

### Covariate file (space/tab delimited)

```
FID IID age sex PC1 PC2 PC3 ...
1001 1001 45 1 0.012 -0.034 0.056 ...
```

### SNP list file (single-column RSIDs)

When using `--snps`, provide a single-column file of RSIDs (one per line), e.g.:

```
rs123
rs456
rs789
```

These RSIDs are matched against the BGEN file's index. Only variants present in both the list and the BGEN are processed.

Note: `chr:pos` formats are not accepted in this file - only RSIDs.

## Output Files

| File | Description |
|------|-------------|
| `{out}.stage1.tsv.gz` | Per-SNP betas and SEs for each tau |
| `{out}.cov.gz` | Binary covariance matrices (upper triangle; int8 by default, float32 optional) |
| `{out}.cov.scale.gz` | Per-SNP float32 scales (only when `--cov-dtype int8`) |
| `{out}.cov.meta.json` | Covariance storage metadata (`cov_dtype`, record layout) |
| `{out}.Rtau.tsv` | T×T correlation matrix (with `--output-rtau`) |

### Stage 1 TSV format

```
snp_id  chr  bp  effect_allele  other_allele  beta_tau_0.10  se_tau_0.10  beta_tau_0.15  se_tau_0.15  ...
rs123   22   1000  A  G  0.0123         0.0045       0.0089         0.0041       ...
```

### Covariance binary format

Default (`--cov-dtype int8`):
- For each SNP, Stage 1 stores the upper triangle of the T×T covariance matrix as int8.
- A per-SNP float32 scale is written to `{out}.cov.scale.gz`.
- Reconstruct as: `cov_upper_float = cov_upper_int8 * scale`.

Legacy mode (`--cov-dtype float32`):
- For each SNP, the upper triangle is stored directly as float32.

For T=17 taus:
- Float32: 153 values = 612 bytes/SNP.
- Int8 + scale: 153 bytes + 4-byte scale = 157 bytes/SNP (~74% smaller).

## Performance

- **C++ kernel**: OpenMP-parallelized block score computation
- **Batched processing**: Configurable batch size for memory efficiency (default 500)
- **BGEN loading**: delay_parsing + index pre-filtering for efficient streaming
- **Speed**: ~800 SNPs/second for full scans, ~22s for 17k HM3 SNPs on chr22

## How It Works

1. **Load phenotype**: Read raw phenotype column from file
2. **Compute RIF**: Build Recentered Influence Function matrix for each tau
3. **Residualize**: Regress out covariates from RIF (not from raw phenotype)
4. **BGEN pre-filter**: Use BGEN index to find valid SNPs from user list
5. **Stream variants**: Iterate BGEN with delay_parsing, extract dosages on demand
6. **Block scores**: Compute per-block sufficient statistics for jackknife
7. **Covariance**: Estimate full T×T covariance per SNP for accurate Stage 2 SEs

## API Reference

### `run_stage1(G, Y, taus, covariates=None, n_blocks=None, seed=42, n_threads=1)`

Full Stage 1 pipeline on arrays.

**Parameters:**
- `G`: (N, M) genotype dosage matrix
- `Y`: (N,) raw phenotype (RIF computed internally)
- `taus`: (T,) quantile levels
- `covariates`: (N, K) optional covariate matrix
- `n_blocks`: number of jackknife blocks
- `seed`: random seed for block assignment
- `n_threads`: OpenMP threads

**Returns:** dict with `betas`, `se`, `covariances`, `q_tau`, `taus`

Important:
- Stage 1 now enforces `n_blocks > n_taus` for covariance stability.
- Recommended for robust Stage 2 `tau_cov`: `n_blocks >= n_taus + 20`.
- A practical default is 45 taus in `[0.02, 0.98]` with auto blocks `max(T + 20, 64)` (i.e., 65 for 45 taus).
