#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 path/to/weights.rds [stage1.tsv.gz] [cov.gz]" >&2
  exit 1
fi

WEIGHTS="$1"
STAGE1="${2:-}"
COV="${3:-}"

Rscript -e '
args <- commandArgs(trailingOnly = TRUE)
weights <- args[1]
stage1 <- if (length(args) >= 2) args[2] else NA
covf <- if (length(args) >= 3) args[3] else NA

options(run_identifiability_examples = FALSE)
source("stage2/diagnose_identifiability.R")

wobj <- readRDS(weights)
if (is.list(wobj) && !is.null(wobj$W)) {
  W <- as.matrix(wobj$W)
  taus <- wobj$taus
} else {
  W <- as.matrix(wobj)
  taus <- NULL
}
if (is.null(taus)) stop("weights must include taus")
if (nrow(W) != length(taus) && ncol(W) == length(taus)) W <- t(W)

if (!is.na(stage1) && !is.na(covf)) {
  emp <- Sigma_tau_from_stage1_cov(stage1, covf, max_snps = 20000)
  Sigma_tau <- emp$Sigma_tau
  taus <- emp$taus
} else {
  Sigma_tau <- make_Sigma_tau(taus, n = 10000)
}

res <- diagnose_identifiability(W, Sigma_tau)
print(res)
' -- "$WEIGHTS" "$STAGE1" "$COV"
