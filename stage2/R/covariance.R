#' Convert upper-triangle covariance to full array
#'
#' @param cov_upper Numeric matrix (n_snps x n_cov_params) of upper-triangle covariances.
#' @param K Integer number of parameters.
#' @param param_names Optional character vector of length K.
#' @return Array of shape (K x K x n_snps).
#' @keywords internal
.param_cov_upper_to_array <- function(cov_upper, K = NULL, param_names = NULL) {
  cov_upper <- as.matrix(cov_upper)
  n_snps <- nrow(cov_upper)
  n_cov <- ncol(cov_upper)
  if (is.null(K)) {
    K <- as.integer((sqrt(8 * n_cov + 1) - 1) / 2)
  }
  if (K * (K + 1) / 2 != n_cov) {
    stop("cov_upper ncol does not match K*(K+1)/2")
  }
  arr <- array(NA_real_, dim = c(K, K, n_snps))
  idx <- 1L
  for (r in seq_len(K)) {
    for (c in r:K) {
      arr[r, c, ] <- cov_upper[, idx]
      arr[c, r, ] <- cov_upper[, idx]
      idx <- idx + 1L
    }
  }
  if (!is.null(param_names)) {
    dimnames(arr) <- list(param_names, param_names, NULL)
  }
  arr
}

#' Derive mashr inputs from parameter covariance
#'
#' @param params Numeric matrix (K x n_snps) of parameter estimates (Bhat).
#' @param cov_upper Numeric matrix (n_snps x n_cov_params) upper-triangle covariances.
#' @param mode One of "array" (per-SNP correlation array) or "mean" (single average correlation matrix).
#' @return List with Bhat, Shat, V.
#' @export
mashr_inputs_from_cov <- function(params, cov_upper, mode = c("array", "mean")) {
  mode <- match.arg(mode)
  params <- as.matrix(params)
  cov_upper <- as.matrix(cov_upper)
  K <- nrow(params)
  if (ncol(cov_upper) != K * (K + 1) / 2) {
    stop("cov_upper has incompatible number of columns for params")
  }
  cov_arr <- .param_cov_upper_to_array(cov_upper, K = K)
  Shat <- matrix(NA_real_, nrow = ncol(params), ncol = K)
  colnames(Shat) <- rownames(params)
  # compute Shat and correlation(s)
  if (mode == "array") {
    V <- array(NA_real_, dim = c(K, K, ncol(params)))
    for (j in seq_len(ncol(params))) {
      cov_j <- cov_arr[, , j]
      sd_j <- sqrt(pmax(diag(cov_j), 0))
      Shat[j, ] <- sd_j
      denom <- outer(sd_j, sd_j)
      V[, , j] <- cov_j / pmax(denom, 1e-12)
    }
  } else {
    # mean correlation across SNPs
    V <- matrix(0, nrow = K, ncol = K)
    for (j in seq_len(ncol(params))) {
      cov_j <- cov_arr[, , j]
      sd_j <- sqrt(pmax(diag(cov_j), 0))
      Shat[j, ] <- sd_j
      denom <- outer(sd_j, sd_j)
      V <- V + cov_j / pmax(denom, 1e-12)
    }
    V <- V / max(ncol(params), 1L)
  }
  list(Bhat = t(params), Shat = Shat, V = V)
}

#' Write parameter covariance upper triangles to a binary file
#'
#' @param cov_upper Numeric matrix (n_snps x n_cov_params).
#' @param path Output path (.cov.gz recommended).
#' @return Invisible TRUE.
#' @keywords internal
.write_param_cov_upper <- function(cov_upper, path) {
  con <- gzfile(path, "wb")
  on.exit(close(con), add = TRUE)
  # write row-by-row in float32
  writeBin(as.numeric(t(cov_upper)), con, size = 4)
  invisible(TRUE)
}

#' Build upper-triangle column names from parameter names
#'
#' @param param_names Character vector.
#' @return Character vector of length K*(K+1)/2.
#' @keywords internal
.param_cov_upper_names <- function(param_names) {
  K <- length(param_names)
  out <- character(K * (K + 1) / 2)
  idx <- 1L
  for (r in seq_len(K)) {
    for (c in r:K) {
      out[idx] <- paste0("cov_", param_names[r], "_", param_names[c])
      idx <- idx + 1L
    }
  }
  out
}
