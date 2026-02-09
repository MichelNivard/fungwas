#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 6 ]]; then
  cat <<USAGE
Usage:
  stage2/scripts/run_full_pipeline_simulation.sh \
    <out_prefix> <weights_arg> <n_samples> <n_snps> <n_threads> <k_eval_csv>

Example:
  stage2/scripts/run_full_pipeline_simulation.sh \
    results/fullsim/bmi_fullsim \
    add=inputs/data/weights/bmi_100tau/weights_additive_bmi.rds,shift=inputs/data/weights/bmi_100tau/weights_multiplicative_bmi.rds,both=inputs/data/weights/bmi_100tau/weights_both_bmi.rds,vqtl=inputs/data/weights/bmi_100tau/weights_vqtl_bmi.rds \
    12000 1500 8 9,13,17,25,33,40,60
USAGE
  exit 1
fi

OUT_PREFIX="$1"
WEIGHTS="$2"
N_SAMPLES="$3"
N_SNPS="$4"
N_THREADS="$5"
K_EVAL="$6"

OUT_DIR="$(dirname "$OUT_PREFIX")"
mkdir -p "$OUT_DIR"

python3 stage2/scripts/simulate_full_pipeline_stage1.py \
  --out-prefix="$OUT_PREFIX" \
  --weights="$WEIGHTS" \
  --n-samples="$N_SAMPLES" \
  --n-snps="$N_SNPS" \
  --n-threads="$N_THREADS"

Rscript stage2/scripts/evaluate_full_pipeline_taucount.R \
  --stage1-tsv="${OUT_PREFIX}.stage1.tsv.gz" \
  --cov-file="${OUT_PREFIX}.cov.gz" \
  --truth-tsv="${OUT_PREFIX}.truth.tsv.gz" \
  --weights="$WEIGHTS" \
  --out-dir="${OUT_PREFIX}.eval" \
  --k-eval="$K_EVAL" \
  --kmax="$(python3 - <<PY
vals=list(map(int,"$K_EVAL".split(',')))
print(max(vals))
PY
)"

echo "Done. Summary: ${OUT_PREFIX}.eval/summary_recovery.tsv"
