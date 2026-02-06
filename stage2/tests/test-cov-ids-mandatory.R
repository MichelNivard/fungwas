# ============================================================
# Test: tau_cov requires covariance IDs sidecar
# ============================================================

library(data.table)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- sub("^--file=", "", file_arg)
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(pkg_root, "R", "gwas.R"))

# Minimal Stage 1 file (2 SNPs, 1 tau)
workdir <- tempfile("cov_ids_test_")
dir.create(workdir)

stage1_file <- file.path(workdir, "test.stage1.tsv")
stage1_dt <- data.table(
  snp_id = c("rs1", "rs2"),
  chr = c(1, 1),
  bp = c(100, 200),
  effect_allele = c("A", "C"),
  other_allele = c("G", "T"),
  beta_tau_0.1 = c(0.0, 0.0),
  se_tau_0.1 = c(1.0, 1.0)
)

fwrite(stage1_dt, stage1_file, sep = "\t")

# Dummy covariance file (sidecar intentionally missing)
cov_file <- file.path(workdir, "test.cov.gz")
con <- gzfile(cov_file, "wb")
writeBin(raw(0), con)
close(con)

# Minimal W (not used before the expected error)
W <- matrix(1, nrow = 1, ncol = 1)

err <- tryCatch({
  param_gwas_from_file(
    stage1_file = stage1_file,
    W = W,
    se_mode = "tau_cov",
    cov_file = cov_file
  )
  NULL
}, error = function(e) e)

if (is.null(err)) {
  stop("Expected error for missing covariance IDs sidecar, but got none")
}

msg <- conditionMessage(err)
if (!grepl("Missing covariance IDs file", msg, fixed = TRUE)) {
  stop("Unexpected error message: ", msg)
}

cat("OK: missing covariance IDs sidecar triggers error\n")
