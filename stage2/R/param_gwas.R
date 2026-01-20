#' Model tournament for parametric GWAS
#'
#' Fits multiple parametric weight matrices to Stage 1 results in one pass
#' and returns per-SNP AIC/BIC values for model selection.
#'
#' @param stage1 Output list from quantile_gwas() containing Q_slope and SE_tau.
#' @param W_list List of weight matrices (each T x K).
#' @param N Sample size used for BIC.
#' @param n_threads Number of OpenMP threads.
#'
#' @return Matrix with AIC/BIC columns per model.
#'
#' @export
param_gwas_multi <- function(stage1, W_list, N, n_threads = 1L) {
  if (is.null(stage1$Q_slope) || is.null(stage1$SE_tau)) {
    stop("stage1 must include Q_slope and SE_tau.")
  }
  if (!is.list(W_list) || length(W_list) == 0) {
    stop("W_list must be a non-empty list of weight matrices.")
  }

  beta_stage1 <- t(stage1$Q_slope)
  se_stage1 <- t(stage1$SE_tau)

  result <- fit_multi_models(beta_stage1, se_stage1, W_list, N, n_threads)

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
