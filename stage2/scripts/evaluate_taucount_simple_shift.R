#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

usage <- function() {
  cat(
"Usage:\n",
"  Rscript stage2/scripts/evaluate_taucount_simple_shift.R \\\n",
"    --stage1-tsv=... --cov-file=... --truth-tsv=... --shift-weight=... --out-tsv=... \\\n",
"    --k-eval=9,13,17,25,33,40,60,80,100\n",
sep = "")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(stage1_tsv=NULL, cov_file=NULL, truth_tsv=NULL, shift_weight=NULL,
              out_tsv="results/simple_shift_k_sweep.tsv", k_eval=NULL,
              tau_min=NA_real_, tau_max=NA_real_)
  for (a in args) {
    if (a %in% c("--help","-h")) { usage(); quit(save="no", status=0) }
    if (!startsWith(a, "--") || !grepl("=", a, fixed=TRUE)) stop("Malformed arg: ", a)
    kv <- strsplit(sub("^--","",a), "=", fixed=TRUE)[[1]]
    key <- kv[1]; val <- paste(kv[-1], collapse="=")
    if (key=="stage1-tsv") opt$stage1_tsv <- val
    else if (key=="cov-file") opt$cov_file <- val
    else if (key=="truth-tsv") opt$truth_tsv <- val
    else if (key=="shift-weight") opt$shift_weight <- val
    else if (key=="out-tsv") opt$out_tsv <- val
    else if (key=="k-eval") opt$k_eval <- as.integer(strsplit(val,",",fixed=TRUE)[[1]])
    else if (key=="tau-min") opt$tau_min <- as.numeric(val)
    else if (key=="tau-max") opt$tau_max <- as.numeric(val)
    else stop("Unknown arg: ", key)
  }
  req <- c("stage1_tsv","cov_file","truth_tsv","shift_weight")
  miss <- req[vapply(req, function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
  if (length(miss)>0) stop("Missing required args: ", paste(miss, collapse=", "))
  opt
}

cov_upper_len <- function(Tt) as.integer(Tt*(Tt+1)/2)

upper_to_sym <- function(upper, Tt) {
  S <- matrix(0, Tt, Tt)
  idx <- 1L
  for (r in seq_len(Tt)) {
    for (c in r:Tt) {
      S[r,c] <- upper[idx]; S[c,r] <- upper[idx]; idx <- idx + 1L
    }
  }
  S
}

extract_cov_rows <- function(cov_file, Tt, row_idx) {
  row_idx <- as.integer(row_idx)
  need <- sort(unique(row_idx))
  n_cov <- cov_upper_len(Tt)
  out_sorted <- matrix(NA_real_, nrow=length(need), ncol=n_cov)
  scale_file <- sub("\\.cov\\.gz$", ".cov.scale.gz", cov_file)
  has_scale <- file.exists(scale_file)
  con <- gzfile(cov_file, "rb")
  on.exit(close(con), add=TRUE)
  con_scale <- NULL
  if (has_scale) {
    con_scale <- gzfile(scale_file, "rb")
    on.exit(close(con_scale), add=TRUE)
  }
  cur <- 1L; k <- 1L
  while (k <= length(need)) {
    if (has_scale) {
      x_q <- readBin(con, "integer", n=n_cov, size=1, signed=TRUE)
      if (length(x_q) < n_cov) stop("int8 cov.gz ended early")
      sc <- readBin(con_scale, "numeric", n=1L, size=4)
      if (length(sc) < 1L) stop("cov.scale.gz ended early")
      x <- as.numeric(x_q) * sc
    } else {
      x <- readBin(con, "numeric", n=n_cov, size=4)
      if (length(x) < n_cov) stop("cov.gz ended early")
    }
    if (cur == need[k]) { out_sorted[k,] <- x; k <- k + 1L }
    cur <- cur + 1L
  }
  out_sorted[match(row_idx, need), , drop=FALSE]
}

fit_gls <- function(b, S, W) {
  S <- (S + t(S))/2
  e <- eigen(S, symmetric=TRUE)
  e$values[e$values < 1e-10] <- 1e-10
  S <- e$vectors %*% (e$values * t(e$vectors))
  C <- tryCatch(chol(S), error=function(e) NULL)
  if (is.null(C)) return(rep(NA_real_, ncol(W)))
  X <- backsolve(C, W, transpose=TRUE)
  z <- backsolve(C, b, transpose=TRUE)
  tryCatch(qr.solve(X, z), error=function(e) rep(NA_real_, ncol(W)))
}

main <- function() {
  opt <- parse_args()
  dir.create(dirname(opt$out_tsv), recursive=TRUE, showWarnings=FALSE)

  s1 <- fread(opt$stage1_tsv)
  truth <- fread(opt$truth_tsv)

  beta_cols <- grep("^beta_tau_", names(s1), value=TRUE)
  taus <- as.numeric(sub("^beta_tau_","", beta_cols))
  ord <- order(taus)
  taus <- taus[ord]
  beta_cols <- beta_cols[ord]

  true_beta_cols <- grep("^true_beta_tau_", names(truth), value=TRUE)
  if (length(true_beta_cols) != length(taus)) stop("truth taus mismatch")

  wobj <- readRDS(opt$shift_weight)
  W <- if (is.list(wobj) && !is.null(wobj$W)) wobj$W else wobj
  wt <- if (is.list(wobj) && !is.null(wobj$taus)) as.numeric(wobj$taus) else NULL
  W <- as.matrix(W)
  if (!is.null(wt)) {
    if (nrow(W) != length(wt) && ncol(W) == length(wt)) W <- t(W)
    idx <- vapply(taus, function(t) { j <- which.min(abs(wt-t)); if (abs(wt[j]-t)>1e-8) NA_integer_ else j }, integer(1))
    if (anyNA(idx)) stop("shift weight taus mismatch")
    W <- W[idx,,drop=FALSE]
  }

  s1[, row_idx := .I]
  dt <- merge(s1[, c("snp_id","row_idx", beta_cols), with=FALSE], truth, by="snp_id", all=FALSE)
  setorder(dt, row_idx)

  cov_rows <- extract_cov_rows(opt$cov_file, length(taus), dt$row_idx)

  if (is.null(opt$k_eval)) opt$k_eval <- c(9L,13L,17L,25L,33L,40L,60L,80L,100L)
  eligible <- seq_along(taus)
  if (is.finite(opt$tau_min)) eligible <- eligible[taus[eligible] >= opt$tau_min]
  if (is.finite(opt$tau_max)) eligible <- eligible[taus[eligible] <= opt$tau_max]
  if (length(eligible) < 2L) stop("Fewer than 2 taus in selected tau range")

  k_eval <- sort(unique(pmin(pmax(opt$k_eval, 2L), length(eligible))))

  out <- vector("list", length(k_eval)); oi <- 1L
  theta2_vals <- vector("list", length(k_eval))

  for (k in k_eval) {
    idx_loc <- unique(round(seq(1, length(eligible), length.out=k)))
    idx_loc <- sort(pmin(pmax(idx_loc, 1L), length(eligible)))
    idx <- eligible[idx_loc]

    mse <- rep(NA_real_, nrow(dt))
    theta_err <- rep(NA_real_, nrow(dt))

    for (i in seq_len(nrow(dt))) {
      b_full <- as.numeric(unlist(dt[i, ..beta_cols]))
      b <- b_full[idx]
      S_full <- upper_to_sym(cov_rows[i,], length(taus))
      S <- S_full[idx, idx, drop=FALSE]
      Wk <- W[idx,,drop=FALSE]

      th <- fit_gls(b, S, Wk)
      if (anyNA(th)) next

      pred_full <- as.numeric(W %*% th)
      true_full <- as.numeric(unlist(dt[i, ..true_beta_cols]))
      mse[i] <- mean((pred_full - true_full)^2)

      # Use the last non-missing true theta column for robust single-model runs.
      if (length(theta_cols <- grep("^theta_", names(dt), value=TRUE)) > 0) {
        tv <- as.numeric(unlist(dt[i, ..theta_cols]))
        keep <- which(is.finite(tv))
        if (length(keep) > 0) {
          p <- keep[length(keep)]
          if (p <= length(th)) theta_err[i] <- th[p] - tv[p]
        }
      }
    }

    out[[oi]] <- data.table(
      k = k,
      tau_min = if (is.finite(opt$tau_min)) opt$tau_min else min(taus),
      tau_max = if (is.finite(opt$tau_max)) opt$tau_max else max(taus),
      n_taus_eligible = length(eligible),
      n = sum(is.finite(mse)),
      median_profile_mse = median(mse, na.rm=TRUE),
      mean_profile_mse = mean(mse, na.rm=TRUE),
      rmse_theta2 = sqrt(mean(theta_err^2, na.rm=TRUE)),
      bias_theta2 = mean(theta_err, na.rm=TRUE)
    )
    oi <- oi + 1L
  }

  res <- rbindlist(out)
  fwrite(res, opt$out_tsv, sep="\t")
  print(res)
}

main()
