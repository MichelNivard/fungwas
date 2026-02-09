#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

usage <- function() {
  cat(
"Usage:\n",
"  Rscript stage2/scripts/evaluate_full_pipeline_taucount.R \\\n",
"    --stage1-tsv path/to/sim.stage1.tsv.gz \\\n",
"    --cov-file path/to/sim.cov.gz \\\n",
"    --truth-tsv path/to/sim.truth.tsv.gz \\\n",
"    --weights add=path/add.rds,shift=path/shift.rds,both=path/both.rds,vqtl=path/vqtl.rds \\\n",
"    --out-dir results/full_pipeline_eval\\n\\n",
"Optional args:\n",
"  --kmax 60\n",
"  --k-eval 9,13,17,25,33,40,60\n",
"  --grid-mode both      # no_cov | emp_cov | both\n",
"  --max-snps 2000\n",
"  --seed 6520\n",
"  --help\n",
sep = "")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(
    stage1_tsv = NULL,
    cov_file = NULL,
    truth_tsv = NULL,
    weights = NULL,
    out_dir = "results/full_pipeline_eval",
    kmax = 60L,
    k_eval = NULL,
    grid_mode = "both",
    max_snps = 2000L,
    seed = 6520L
  )

  for (a in args) {
    if (a %in% c("--help", "-h")) {
      usage()
      quit(save = "no", status = 0)
    }
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) {
      stop("Malformed argument: ", a)
    }
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- paste(kv[-1], collapse = "=")

    if (key == "stage1-tsv") opt$stage1_tsv <- val
    else if (key == "cov-file") opt$cov_file <- val
    else if (key == "truth-tsv") opt$truth_tsv <- val
    else if (key == "weights") opt$weights <- val
    else if (key == "out-dir") opt$out_dir <- val
    else if (key == "kmax") opt$kmax <- as.integer(val)
    else if (key == "k-eval") opt$k_eval <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key == "grid-mode") opt$grid_mode <- trimws(tolower(val))
    else if (key == "max-snps") opt$max_snps <- as.integer(val)
    else if (key == "seed") opt$seed <- as.integer(val)
    else stop("Unknown argument: --", key)
  }

  req <- c("stage1_tsv", "cov_file", "truth_tsv", "weights")
  miss <- req[vapply(req, function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
  if (length(miss) > 0) stop("Missing required args: ", paste(miss, collapse = ", "))

  if (!opt$grid_mode %in% c("no_cov", "emp_cov", "both")) {
    stop("--grid-mode must be one of: no_cov, emp_cov, both")
  }

  opt
}

parse_weights_arg <- function(x) {
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  out <- vector("list", length(parts))
  names(out) <- rep("", length(parts))
  for (i in seq_along(parts)) {
    kv <- strsplit(parts[i], "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Bad weights entry: ", parts[i])
    out[[i]] <- trimws(kv[2])
    names(out)[i] <- trimws(kv[1])
  }
  out
}

read_weight_aligned <- function(path, taus_ref) {
  w <- readRDS(path)
  W <- if (is.list(w) && !is.null(w$W)) w$W else w
  wt <- if (is.list(w) && !is.null(w$taus)) as.numeric(w$taus) else NULL
  W <- as.matrix(W)

  if (!is.null(wt)) {
    if (nrow(W) != length(wt) && ncol(W) == length(wt)) W <- t(W)
    idx <- vapply(taus_ref, function(tau) {
      d <- abs(wt - tau)
      j <- which.min(d)
      if (d[j] > 1e-8) NA_integer_ else j
    }, integer(1))
    if (anyNA(idx)) stop("Weight taus mismatch: ", path)
    W <- W[idx, , drop = FALSE]
  } else {
    if (nrow(W) == length(taus_ref)) {
      # pass
    } else if (ncol(W) == length(taus_ref)) {
      W <- t(W)
    } else {
      stop("Weight dimensions incompatible: ", path)
    }
  }

  storage.mode(W) <- "double"
  W
}

proj <- function(M) {
  M <- as.matrix(M)
  M %*% solve(crossprod(M), t(M))
}

build_contrast_rows <- function(W_list) {
  P <- lapply(W_list, proj)
  mods <- names(P)
  out <- list()
  k <- 1L
  for (i in seq_len(length(mods) - 1L)) {
    for (j in (i + 1L):length(mods)) {
      out[[k]] <- P[[mods[i]]] - P[[mods[j]]]
      names(out)[k] <- paste0(mods[i], "_vs_", mods[j])
      k <- k + 1L
    }
  }
  out
}

build_contrast_rows_empcov <- function(W_list, R_tau) {
  U <- chol(R_tau)
  U_inv <- backsolve(U, diag(ncol(U)))
  W_pw <- lapply(W_list, function(W) U_inv %*% W)
  build_contrast_rows(W_pw)
}

obj_fixed_monotone <- function(S, D_list) {
  if (length(S) == 0L) return(0)
  val <- 0
  for (D in D_list) {
    X <- D[S, , drop = FALSE]
    G <- crossprod(X)
    val <- val + as.numeric(determinant(diag(nrow(G)) + G, logarithm = TRUE)$modulus[1])
  }
  val
}

run_greedy <- function(taus, D_list, Kmax, method) {
  Kmax <- min(as.integer(Kmax), length(taus))
  selected <- integer(0)
  cur <- 0
  rows <- vector("list", Kmax)
  for (step in seq_len(Kmax)) {
    cand <- setdiff(seq_along(taus), selected)
    gains <- vapply(cand, function(ci) {
      obj_fixed_monotone(c(selected, ci), D_list) - cur
    }, numeric(1))
    j <- which.max(gains)
    pick <- cand[j]
    selected <- c(selected, pick)
    cur <- obj_fixed_monotone(selected, D_list)
    rows[[step]] <- data.table(
      step = step,
      add_idx = pick,
      add_tau = taus[pick],
      marginal_gain = gains[j],
      cum_objective = cur,
      method = method
    )
  }
  out <- rbindlist(rows)
  out[, rel_gain := marginal_gain / marginal_gain[1]]
  out[, cum_frac := cum_objective / max(cum_objective)]
  out
}

cov_upper_len <- function(Tt) as.integer(Tt * (Tt + 1) / 2)

upper_to_sym <- function(upper, Tt) {
  S <- matrix(0, Tt, Tt)
  idx <- 1L
  for (r in seq_len(Tt)) {
    for (c in r:Tt) {
      S[r, c] <- upper[idx]
      S[c, r] <- upper[idx]
      idx <- idx + 1L
    }
  }
  S
}

extract_cov_rows <- function(cov_file, Tt, row_idx) {
  row_idx <- as.integer(row_idx)
  need <- sort(unique(row_idx))
  n_cov <- cov_upper_len(Tt)
  out_sorted <- matrix(NA_real_, nrow = length(need), ncol = n_cov)
  scale_file <- sub("\\.cov\\.gz$", ".cov.scale.gz", cov_file)
  has_scale <- file.exists(scale_file)

  con <- gzfile(cov_file, "rb")
  on.exit(close(con), add = TRUE)
  con_scale <- NULL
  if (has_scale) {
    con_scale <- gzfile(scale_file, "rb")
    on.exit(close(con_scale), add = TRUE)
  }

  cur <- 1L
  k <- 1L
  while (k <= length(need)) {
    if (has_scale) {
      x_q <- readBin(con, "integer", n = n_cov, size = 1, signed = TRUE)
      if (length(x_q) < n_cov) stop("Reached end of int8 cov file early")
      sc <- readBin(con_scale, "numeric", n = 1L, size = 4)
      if (length(sc) < 1L) stop("Reached end of cov scale file early")
      x <- as.numeric(x_q) * sc
    } else {
      x <- readBin(con, "numeric", n = n_cov, size = 4)
      if (length(x) < n_cov) stop("Reached end of cov file early")
    }
    if (cur == need[k]) {
      out_sorted[k, ] <- x
      k <- k + 1L
    }
    cur <- cur + 1L
  }

  out_sorted[match(row_idx, need), , drop = FALSE]
}

fit_dwls <- function(b, se, W) {
  se <- pmax(se, 1e-12)
  X <- W * (1 / se)
  z <- b / se
  theta <- tryCatch(qr.solve(X, z), error = function(e) rep(NA_real_, ncol(W)))
  if (anyNA(theta)) return(list(theta = theta, Q = NA_real_, AIC = NA_real_))
  r <- b - as.numeric(W %*% theta)
  Q <- sum((r / se)^2)
  list(theta = theta, Q = Q, AIC = Q + 2 * ncol(W))
}

fit_gls <- function(b, S, W) {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  e$values[e$values < 1e-10] <- 1e-10
  S <- e$vectors %*% (e$values * t(e$vectors))

  C <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(C)) return(list(theta = rep(NA_real_, ncol(W)), Q = NA_real_, AIC = NA_real_))

  X <- backsolve(C, W, transpose = TRUE)
  z <- backsolve(C, b, transpose = TRUE)
  theta <- tryCatch(qr.solve(X, z), error = function(e) rep(NA_real_, ncol(W)))
  if (anyNA(theta)) return(list(theta = theta, Q = NA_real_, AIC = NA_real_))

  r <- b - as.numeric(W %*% theta)
  Q <- tryCatch(as.numeric(crossprod(r, solve(S, r))), error = function(e) NA_real_)
  list(theta = theta, Q = Q, AIC = Q + 2 * ncol(W))
}

main <- function() {
  opt <- parse_args()
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$seed)

  s1 <- fread(opt$stage1_tsv)
  truth <- fread(opt$truth_tsv)

  beta_cols <- grep("^beta_tau_", names(s1), value = TRUE)
  se_cols <- grep("^se_tau_", names(s1), value = TRUE)
  taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
  ord <- order(taus)
  taus <- taus[ord]
  beta_cols <- beta_cols[ord]
  se_cols <- se_cols[ord]

  s1[, row_idx := .I]
  dt <- merge(
    s1[, c("snp_id", "row_idx", beta_cols, se_cols), with = FALSE],
    truth,
    by = "snp_id",
    all = FALSE
  )
  if (nrow(dt) == 0) stop("No overlapping SNPs between stage1 and truth")

  if (!is.na(opt$max_snps) && opt$max_snps > 0 && nrow(dt) > opt$max_snps) {
    dt <- dt[sample(.N, opt$max_snps)]
  }
  setorder(dt, row_idx)

  cov_rows <- extract_cov_rows(opt$cov_file, length(taus), dt$row_idx)

  w_map <- parse_weights_arg(opt$weights)
  W_list <- lapply(w_map, read_weight_aligned, taus_ref = taus)
  model_names <- names(W_list)

  D_nocov <- build_contrast_rows(W_list)
  grids <- list()
  if (opt$grid_mode %in% c("no_cov", "both")) {
    grids$no_cov <- run_greedy(taus, D_nocov, opt$kmax, "no_cov")
  }
  if (opt$grid_mode %in% c("emp_cov", "both")) {
    B <- as.matrix(dt[, ..beta_cols])
    R_emp <- cor(B, use = "pairwise.complete.obs")
    e <- eigen(R_emp, symmetric = TRUE)
    e$values[e$values < 1e-10] <- 1e-10
    R_emp <- e$vectors %*% (e$values * t(e$vectors))
    D_empcov <- build_contrast_rows_empcov(W_list, R_emp)
    grids$emp_cov <- run_greedy(taus, D_empcov, opt$kmax, "emp_cov")
  }

  if (length(grids) == 0) stop("No grids computed")

  k_eval <- opt$k_eval
  if (is.null(k_eval) || length(k_eval) == 0L) k_eval <- sort(unique(c(9L, 13L, 17L, 25L, 33L, 40L, opt$kmax)))

  true_beta_cols <- grep("^true_beta_tau_", names(dt), value = TRUE)
  if (length(true_beta_cols) != length(taus)) stop("Truth file must have true_beta_tau_* for all taus")

  theta_cols <- grep("^theta_", names(dt), value = TRUE)

  details <- list(); di <- 1L
  theta_rows <- list(); ti <- 1L

  for (grid_name in names(grids)) {
    gtab <- grids[[grid_name]]
    kk <- sort(unique(pmin(pmax(k_eval, 1L), nrow(gtab))))

    for (k in kk) {
      idx <- gtab$add_idx[seq_len(k)]

      for (se_mode in c("dwls", "tau_cov")) {
        for (i in seq_len(nrow(dt))) {
          b_full <- as.numeric(unlist(dt[i, ..beta_cols]))
          se_full <- as.numeric(unlist(dt[i, ..se_cols]))
          b <- b_full[idx]
          se <- se_full[idx]

          S_full <- upper_to_sym(cov_rows[i, ], length(taus))
          S <- S_full[idx, idx, drop = FALSE]

          fits <- lapply(model_names, function(m) {
            W <- W_list[[m]][idx, , drop = FALSE]
            if (se_mode == "dwls") fit_dwls(b, se, W) else fit_gls(b, S, W)
          })
          names(fits) <- model_names

          aic <- vapply(fits, function(x) x$AIC, numeric(1))
          ordm <- order(aic)
          winner <- model_names[ordm[1]]

          true_model <- dt$true_model[i]
          fit_true <- fits[[true_model]]

          W_win_full <- W_list[[winner]]
          W_true_full <- W_list[[true_model]]

          pred_win_full <- as.numeric(W_win_full %*% fits[[winner]]$theta)
          pred_true_full <- as.numeric(W_true_full %*% fit_true$theta)
          beta_true_full <- as.numeric(unlist(dt[i, ..true_beta_cols]))

          mse_win <- mean((pred_win_full - beta_true_full)^2)
          mse_true <- mean((pred_true_full - beta_true_full)^2)

          details[[di]] <- data.table(
            snp_id = dt$snp_id[i],
            grid = grid_name,
            k = k,
            se_mode = se_mode,
            true_model = true_model,
            winner = winner,
            correct = as.integer(winner == true_model),
            delta = if (length(ordm) >= 2) aic[ordm[2]] - aic[ordm[1]] else NA_real_,
            mse_profile_winner = mse_win,
            mse_profile_true_model_fit = mse_true
          )
          di <- di + 1L

          theta_true <- as.numeric(unlist(dt[i, ..theta_cols]))
          theta_hat <- fit_true$theta
          L <- min(length(theta_hat), length(theta_true))
          if (L > 0) {
            for (p in seq_len(L)) {
              if (is.finite(theta_true[p])) {
                theta_rows[[ti]] <- data.table(
                  snp_id = dt$snp_id[i],
                  grid = grid_name,
                  k = k,
                  se_mode = se_mode,
                  true_model = true_model,
                  param_idx = p,
                  theta_true = theta_true[p],
                  theta_hat = theta_hat[p],
                  theta_err = theta_hat[p] - theta_true[p]
                )
                ti <- ti + 1L
              }
            }
          }
        }
      }
    }
  }

  detail <- rbindlist(details, fill = TRUE)
  theta_dt <- if (length(theta_rows) > 0) rbindlist(theta_rows, fill = TRUE) else data.table()

  sum_acc <- detail[, .(
    n = .N,
    accuracy = mean(correct, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    median_mse_winner = median(mse_profile_winner, na.rm = TRUE),
    median_mse_truefit = median(mse_profile_true_model_fit, na.rm = TRUE)
  ), by = .(grid, k, se_mode)][order(grid, k, se_mode)]

  confusion <- detail[, .N, by = .(grid, k, se_mode, true_model, winner)][order(grid, k, se_mode, true_model, -N)]

  if (nrow(theta_dt) > 0) {
    theta_summary <- theta_dt[, .(
      n = .N,
      bias = mean(theta_err, na.rm = TRUE),
      mse = mean(theta_err^2, na.rm = TRUE),
      rmse = sqrt(mean(theta_err^2, na.rm = TRUE)),
      cor_true_hat = suppressWarnings(cor(theta_true, theta_hat, use = "pairwise.complete.obs"))
    ), by = .(grid, k, se_mode, true_model, param_idx)][order(grid, k, se_mode, true_model, param_idx)]
  } else {
    theta_summary <- data.table()
  }

  fwrite(rbindlist(grids, idcol = "grid", use.names = TRUE), file.path(opt$out_dir, "grid_sequences.tsv"), sep = "\t")
  fwrite(detail, file.path(opt$out_dir, "detail_recovery.tsv.gz"), sep = "\t")
  fwrite(sum_acc, file.path(opt$out_dir, "summary_recovery.tsv"), sep = "\t")
  fwrite(confusion, file.path(opt$out_dir, "confusion_counts.tsv"), sep = "\t")
  fwrite(theta_summary, file.path(opt$out_dir, "theta_bias_mse.tsv"), sep = "\t")

  p1 <- ggplot(sum_acc, aes(x = k, y = accuracy, color = se_mode)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.4) + facet_wrap(~grid) +
    labs(title = "Model recovery accuracy vs tau count", x = "k", y = "Accuracy", color = "SE mode") +
    theme_minimal(base_size = 12)
  ggsave(file.path(opt$out_dir, "accuracy_vs_k.png"), p1, width = 9, height = 5, dpi = 170)

  p2 <- ggplot(sum_acc, aes(x = k, y = median_mse_winner, color = se_mode)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.3) + facet_wrap(~grid) +
    labs(title = "Winner profile MSE vs tau count", x = "k", y = "Median MSE", color = "SE mode") +
    theme_minimal(base_size = 12)
  ggsave(file.path(opt$out_dir, "winner_mse_vs_k.png"), p2, width = 9, height = 5, dpi = 170)

  message("Saved outputs in: ", opt$out_dir)
  print(sum_acc)
}

main()
