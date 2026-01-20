# FungWas Stage 1: Fast RIF Quantile GWAS

Python/C++ implementation for blazing-fast RIF quantile regression with jackknife covariance estimation.

## Installation

```bash
# With conda (recommended)
conda env create -f environment.yml
conda activate fungwas-stage1
pip install .

# Or pip only
pip install .
```

### Conda environment notes (BGEN support)

Stage 1 relies on the Python `bgen` module for BGEN streaming, so the conda
environment installs it via pip (currently pinned to `bgen==1.9.7`). This avoids
HPC module requirements and keeps the CLI self-contained.

```bash
conda env update -n fungwas-stage1 -f environment.yml
```

**Requirement:** `bgenix` must be available for SNP extraction. Install via:
```bash
# On HPC, typically available as module
module load bgenix

# Or download from https://enkre.net/cgi-bin/code/bgen/dir?ci=tip
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
taus = np.arange(0.10, 0.91, 0.05)

# Run Stage 1 - RIF is computed internally
result = run_stage1(G, Y, taus, n_blocks=25, n_threads=4)

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
    --taus "0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90" \
    --blocks 25 \
    --batch-size 1000 \
    --threads 8 \
    --bgenix-path /path/to/bgenix \
    --output-rtau \
    --out results/chr22
```

#### Multi-phenotype CLI (single BGEN pass)

Provide comma-separated phenotype columns to scan multiple phenotypes with one
genotype pass. Outputs are suffixed with the phenotype name.

```bash
fungwas-stage1 \
    --bgen data.chr22.bgen \
    --pheno phenotypes.txt \
    --pheno-cols "phenotype_a,phenotype_b" \
    --covar covariates.txt \
    --snps chr22_hm3.txt \
    --taus "0.10,0.50,0.90" \
    --blocks 25 \
    --batch-size 1000 \
    --threads 8 \
    --bgenix-path /path/to/bgenix \
    --out results/chr22
```

If `--pheno-cols` contains a single phenotype, Stage 1 automatically falls back
to the single-phenotype kernel to avoid extra overhead.

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

### Phenotype file (space-delimited)

Standard format with FID, IID, and raw phenotype column(s):
```
FID IID testosterone age_at_assessment
1001 1001 15.3 45
1002 1002 12.1 52
```

The `--pheno-col` argument specifies which column to analyze. FungWas computes RIF internally.

### Covariate file (space-delimited)

```
FID IID age sex PC1 PC2 PC3 ...
1001 1001 45 1 0.012 -0.034 0.056 ...
```

## Output Files

| File | Description |
|------|-------------|
| `{out}.stage1.tsv.gz` | Per-SNP betas and SEs for each tau |
| `{out}.cov.gz` | Binary covariance matrices (upper triangle, float32) |
| `{out}.Rtau.tsv` | T×T correlation matrix (with `--output-rtau`) |

### Stage 1 TSV format

```
snp_id  chr  bp  a1  a2  beta_tau_0.10  se_tau_0.10  beta_tau_0.15  se_tau_0.15  ...
rs123   22   1000  A  G  0.0123         0.0045       0.0089         0.0041       ...
```

### Covariance binary format

For each SNP, the upper triangle of the T×T covariance matrix is stored as float32.
For T=17 taus: 153 floats = 612 bytes per SNP.

## Performance

- **C++ kernel**: OpenMP-parallelized block score computation
- **Batched processing**: Configurable batch size for memory efficiency
- **Speed**: ~1000 SNPs/second on single thread, scales with `--threads`

## How It Works

1. **Load phenotype**: Read raw phenotype column from file
2. **Compute RIF**: Build Recentered Influence Function matrix for each tau
3. **Residualize**: Regress out covariates from RIF (not from raw phenotype)
4. **Block scores**: Compute per-block sufficient statistics for jackknife
5. **Covariance**: Estimate full T×T covariance per SNP for accurate Stage 2 SEs

## API Reference

### `run_stage1(G, Y, taus, covariates=None, n_blocks=25, seed=42, n_threads=1)`

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
