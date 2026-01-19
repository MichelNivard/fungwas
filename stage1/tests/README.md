# Stage 1 Validation Chain

This folder contains sequential proofs that the jackknife standard errors
match analytic OLS results under progressively more complex models.

## Order of checks

1) `test_jk_ols_equivalence.py`
   - Standard OLS: `Y ~ G` using block jackknife vs analytic OLS.
   - Uses the same preprocessing as the C++ kernel (imputation, MAC/min-var,
     and covariate residualization via QR).

2) `test_jk_equivalence.py`
   - Stage 1 RIF/quantile GWAS: jackknife SEs vs analytic RIF OLS.
   - Uses the same preprocessing as the C++ kernel.

3) Stage 2 comparison (in the R package)
   - See `stage2/tests/proof_stage2_se.R` for Stage 2 SE mode comparisons.

## Running

From the repo root:

```bash
PYTHONPATH=stage1 python stage1/tests/test_jk_ols_equivalence.py
PYTHONPATH=stage1 python stage1/tests/test_jk_equivalence.py
Rscript stage2/tests/proof_stage2_se.R
```
