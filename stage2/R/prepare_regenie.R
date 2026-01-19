#' Prepare REGENIE inputs for RIF regression
#'
#' Computes RIF-transformed phenotypes for a set of taus and writes
#' phenotype and covariate files compatible with REGENIE step 2.
#'
#' @param pheno_df Data frame with FID/IID columns and phenotype columns.
#' @param covar_df Data frame with FID/IID columns and covariate columns.
#' @param phenotypes Character vector of phenotype column names.
#' @param taus Numeric vector of quantile levels.
#' @param covar_cols Character vector of covariate column names to include.
#' @param out_inputs Output directory for REGENIE files.
#' @param out_rtau Optional output directory for R_tau matrices.
#' @param density_floor Minimum density value for RIF computation.
#'
#' @return Invisibly returns a list of outputs for each phenotype.
#' @export
prepare_regenie <- function(
  pheno_df,
  covar_df,
  phenotypes,
  taus,
  covar_cols,
  out_inputs,
  out_rtau = NULL,
  density_floor = 1e-8
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required for prepare_regenie")
  }

  pheno_dt <- data.table::as.data.table(pheno_df)
  covar_dt <- data.table::as.data.table(covar_df)

  if (!all(c("FID", "IID") %in% names(pheno_dt))) {
    stop("pheno_df must include FID and IID columns")
  }
  if (!all(c("FID", "IID") %in% names(covar_dt))) {
    stop("covar_df must include FID and IID columns")
  }
  if (!all(phenotypes %in% names(pheno_dt))) {
    stop("Missing phenotype columns in pheno_df")
  }
  if (!all(covar_cols %in% names(covar_dt))) {
    stop("Missing covariate columns in covar_df")
  }

  data.table::setkeyv(pheno_dt, c("FID", "IID"))
  data.table::setkeyv(covar_dt, c("FID", "IID"))

  out_inputs <- normalizePath(out_inputs, mustWork = FALSE)
  dir.create(out_inputs, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(out_rtau)) {
    out_rtau <- normalizePath(out_rtau, mustWork = FALSE)
    dir.create(out_rtau, recursive = TRUE, showWarnings = FALSE)
  }

  outputs <- list()

  for (pheno_name in phenotypes) {
    merged_df <- merge(
      as.data.frame(pheno_dt),
      as.data.frame(covar_dt),
      by = c("FID", "IID")
    )
    cols_needed <- c("FID", "IID", pheno_name, covar_cols)
    cols_needed <- as.character(unlist(cols_needed))
    merged_df <- merged_df[stats::complete.cases(merged_df), cols_needed, drop = FALSE]
    merged <- data.table::as.data.table(merged_df)

    y <- merged[[pheno_name]]
    rif <- .build_rif_matrix(y, taus, density_floor = density_floor)
    rif_mat <- rif$R
    colnames(rif_mat) <- gsub("\\.", "_", sprintf("RIF_tau%.2f", taus))

    pheno_out <- data.table::data.table(FID = merged$FID, IID = merged$IID)
    pheno_out <- cbind(pheno_out, rif_mat)

    covar_out <- merged[, c("FID", "IID", covar_cols), with = FALSE]

    pheno_path <- file.path(out_inputs, sprintf("regenie_pheno_%s.txt", pheno_name))
    covar_path <- file.path(out_inputs, "regenie_covar.txt")

    data.table::fwrite(pheno_out, pheno_path, sep = " ", na = "NA", quote = FALSE)
    data.table::fwrite(covar_out, covar_path, sep = " ", na = "NA", quote = FALSE)

    rtau_path <- NULL
    if (!is.null(out_rtau)) {
      covar_mat <- as.matrix(merged[, covar_cols, with = FALSE])
      covar_mat <- cbind(1, covar_mat)
      q_mat <- qr.Q(qr(covar_mat))
      rif_resid <- rif_mat - q_mat %*% crossprod(q_mat, rif_mat)
      r_tau <- stats::cor(rif_resid)

      if (requireNamespace("Matrix", quietly = TRUE)) {
        r_tau <- as.matrix(Matrix::nearPD(r_tau, corr = TRUE)$mat)
      }

      rtau_path <- file.path(out_rtau, sprintf("R_tau_%s.tsv", pheno_name))
      data.table::fwrite(as.data.table(r_tau), rtau_path, sep = "\t", quote = FALSE)
    }

    outputs[[pheno_name]] <- list(
      pheno_file = pheno_path,
      covar_file = covar_path,
      rtau_file = rtau_path
    )
  }

  invisible(outputs)
}
