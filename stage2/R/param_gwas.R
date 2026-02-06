#' Model tournament for parametric GWAS
#'
#' Fits multiple parametric weight matrices to Stage 1 results in one pass
#' and returns per-SNP AIC/BIC values for model selection.
#'
#' @param stage1 Output list from quantile_gwas() containing Q_slope and SE_tau.
#' @param W_list List of weight matrices (each T x K).
#' @param N Sample size used for BIC.
#' @param n_threads Number of OpenMP threads.
#' @param se_mode One of \code{"tau_cov"} or \code{"dwls"} (or \code{"diagonal"}).
#'   \code{"tau_cov"} uses per-SNP covariance (recommended when available).
#' @param cov_data Optional list with \code{cov_vec} and \code{offsets} for
#'   \code{se_mode = "tau_cov"}.
#'
#' @return Matrix with AIC/BIC columns per model.
#'
#' @export
param_gwas_multi <- function(stage1, W_list, N, n_threads = 1L,
                             se_mode = c("tau_cov", "dwls", "diagonal"),
                             cov_data = NULL) {
  se_mode <- match.arg(se_mode)
  if (is.null(stage1$Q_slope) || is.null(stage1$SE_tau)) {
    stop("stage1 must include Q_slope and SE_tau.")
  }
  if (!is.list(W_list) || length(W_list) == 0) {
    stop("W_list must be a non-empty list of weight matrices.")
  }

  beta_stage1 <- t(stage1$Q_slope)
  se_stage1 <- t(stage1$SE_tau)

  if (se_mode == "tau_cov") {
    if (is.null(cov_data) || is.null(cov_data$cov_vec) || is.null(cov_data$offsets)) {
      stop("se_mode = 'tau_cov' requires cov_data with cov_vec and offsets.")
    }
    result <- fit_multi_models_tau_cov(
      beta_stage1, cov_data$cov_vec, cov_data$offsets, W_list, N, n_threads
    )
  } else {
    result <- fit_multi_models(beta_stage1, se_stage1, W_list, N, n_threads)
  }

  model_names <- names(W_list)
  if (is.null(model_names) || any(model_names == "")) {
    model_names <- paste0("model", seq_along(W_list))
  }
  colnames(result) <- as.vector(rbind(
    paste0("AIC_", model_names),
    paste0("BIC_", model_names)
  ))

  result
}
