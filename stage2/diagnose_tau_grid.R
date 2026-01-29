# Tau-grid diagnostic for Stage 2 identifiability
#
# Usage (example):
#   Rscript stage2/diagnose_tau_grid.R \
#     --weights inputs/data/weights/weights_both_bmi.rds \
#     --stage1 results/fungwas/combined/bmi.stage1.tsv.gz \
#     --cov results/fungwas/combined/bmi.cov.gz \
#     --grid 9,13,17,25 \
#     --max-snps 20000

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
})

options(run_identifiability_examples = FALSE)
source("stage2/diagnose_identifiability.R")

read_weight <- function(path) {
  obj <- readRDS(path)
  if (is.list(obj) && !is.null(obj$W)) {
    W <- as.matrix(obj$W)
    taus <- obj$taus
  } else {
    W <- as.matrix(obj)
    taus <- NULL
  }
  list(W = W, taus = taus)
}

subset_grid <- function(taus, n) {
  if (n >= length(taus)) return(seq_along(taus))
  idx <- round(seq(1, length(taus), length.out = n))
  unique(pmin(pmax(idx, 1), length(taus)))
}

summarize_grid <- function(W, taus, Sigma_tau) {
  diag <- diagnose_identifiability(W, Sigma_tau)
  list(
    condition_number = diag$condition_number,
    mean_se_inflation = diag$mean_se_inflation,
    max_se_reduction = diag$max_se_reduction
  )
}

option_list <- list(
  make_option(c("--weights"), type = "character", help = "Path to weight RDS (W and taus)", metavar = "FILE"),
  make_option(c("--stage1"), type = "character", default = NULL, help = "Stage1 TSV (for empirical taus)"),
  make_option(c("--cov"), type = "character", default = NULL, help = "Stage1 cov.gz (for empirical Sigma_tau)"),
  make_option(c("--grid"), type = "character", default = "9,13,17,25", help = "Comma-separated grid sizes"),
  make_option(c("--max-snps"), type = "integer", default = 20000, help = "Max SNPs to average for empirical Sigma_tau"),
  make_option(c("--n"), type = "integer", default = 10000, help = "Sample size for theoretical Sigma_tau")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$weights)) stop("--weights is required")

wobj <- read_weight(opt$weights)
if (is.null(wobj$taus)) stop("Weight file must include taus")

if (!is.null(opt$stage1) && !is.null(opt$cov)) {
  emp <- Sigma_tau_from_stage1_cov(opt$stage1, opt$cov, max_snps = opt$`max-snps`)
  Sigma_tau_full <- emp$Sigma_tau
  taus_full <- emp$taus
  if (!all.equal(taus_full, wobj$taus, tolerance = 1e-8)) {
    stop("Stage1 taus do not match weight taus")
  }
} else {
  taus_full <- wobj$taus
  Sigma_tau_full <- make_Sigma_tau(taus_full, n = opt$n)
}

W_full <- wobj$W
if (nrow(W_full) != length(taus_full) && ncol(W_full) == length(taus_full)) {
  W_full <- t(W_full)
}

sizes <- as.integer(strsplit(opt$grid, ",")[[1]])

cat("tau_grid\tcondition_number\tmean_se_inflation\tmax_se_reduction\n")
for (n in sizes) {
  idx <- subset_grid(taus_full, n)
  W_sub <- W_full[idx, , drop = FALSE]
  Sigma_sub <- Sigma_tau_full[idx, idx, drop = FALSE]
  s <- summarize_grid(W_sub, taus_full[idx], Sigma_sub)
  cat(paste(n, sprintf("%.2f", s$condition_number), sprintf("%.2f", s$mean_se_inflation),
            sprintf("%.2f", s$max_se_reduction), sep = "\t"), "\n")
}
