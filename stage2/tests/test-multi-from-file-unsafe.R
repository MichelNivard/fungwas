# ============================================================
# Test: param_gwas_multi_from_file requires cov IDs unless unsafe
# ============================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rcpp))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
this_file <- sub("^--file=", "", file_arg)
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))
source(file.path(pkg_root, "R", "gwas.R"))
source(file.path(pkg_root, "R", "param_gwas.R"))
source(file.path(pkg_root, "R", "RcppExports.R"))
Rcpp::sourceCpp(file.path(pkg_root, "src", "models.cpp"))

workdir <- tempfile("multi_from_file_")
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

cov_file <- file.path(workdir, "test.cov.gz")
con <- gzfile(cov_file, "wb")
writeBin(raw(0), con)
close(con)

W_list <- list(add = matrix(1, nrow = 1, ncol = 1))

err <- tryCatch({
  param_gwas_multi_from_file(
    stage1_file = stage1_file,
    W_list = W_list,
    N = 100,
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

# Now proceed with unsafe flag (should warn but not error)
suppressWarnings(
  param_gwas_multi_from_file(
    stage1_file = stage1_file,
    W_list = W_list,
    N = 100,
    se_mode = "tau_cov",
    cov_file = cov_file,
    unsafe_skip_cov_ids = TRUE
  )
)

cat("OK: unsafe_skip_cov_ids bypasses missing cov IDs check\n")
