#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

usage <- function() {
  cat(
"Usage:\n",
"  Rscript stage2/scripts/taugrid_model_compare_fixed_nocov.R \\\n",
"    --stage1-tsv path/to/stage1.tsv.gz \\\n",
"    --weights add=path/add.rds,shift=path/shift.rds,both=path/both.rds,vqtl=path/vqtl.rds \\\n",
"    --out-dir results/taugrid_generic\\n\\n",
"Alternative input:\n",
"  --stage1-rds path/to/stage1.rds\n",
"    Expected list elements: taus, Q_slope (T x P), SE_tau (T x P).\n\\n",
"Optional args:\n",
"  --kmax 40\n",
"  --k-eval 9,13,17,25,33,40\n",
"  --greedy-mode both    # no_cov | emp_cov | both\n",
"  --max-snps 2000\n",
"  --seed 6520\n",
"  --help\n",
sep = "")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    stage1_tsv = NULL,
    stage1_rds = NULL,
    weights = NULL,
    out_dir = "results/taugrid_generic",
    kmax = 40L,
    k_eval = NULL,
    greedy_mode = "both",
    max_snps = 2000L,
    seed = 6520L
  )

  if (length(args) == 0L) return(out)
  for (a in args) {
    if (a %in% c("--help", "-h")) {
      usage()
      quit(save = "no", status = 0)
    }
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) {
      stop("Malformed argument: ", a, " (expected --key=value)")
    }
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- paste(kv[-1], collapse = "=")

    if (key == "stage1-tsv") out$stage1_tsv <- val
    else if (key == "stage1-rds") out$stage1_rds <- val
    else if (key == "weights") out$weights <- val
    else if (key == "out-dir") out$out_dir <- val
    else if (key == "kmax") out$kmax <- as.integer(val)
    else if (key == "k-eval") out$k_eval <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key == "greedy-mode") out$greedy_mode <- trimws(tolower(val))
    else if (key == "max-snps") out$max_snps <- as.integer(val)
    else if (key == "seed") out$seed <- as.integer(val)
    else stop("Unknown argument: --", key)
  }

  out
}

read_stage1 <- function(stage1_tsv = NULL, stage1_rds = NULL) {
  if (!is.null(stage1_tsv) && !is.null(stage1_rds)) {
    stop("Provide only one of --stage1-tsv or --stage1-rds")
  }
  if (is.null(stage1_tsv) && is.null(stage1_rds)) {
    stop("Provide one of --stage1-tsv or --stage1-rds")
  }

  if (!is.null(stage1_tsv)) {
    s1 <- fread(stage1_tsv)
    beta_cols <- grep("^beta_tau_", names(s1), value = TRUE)
    se_cols <- grep("^se_tau_", names(s1), value = TRUE)
    if (length(beta_cols) == 0 || length(se_cols) == 0) {
      stop("stage1-tsv must include beta_tau_* and se_tau_* columns")
    }

    taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
    ord <- order(taus)
    taus <- taus[ord]
    beta_cols <- beta_cols[ord]
    se_cols <- se_cols[ord]

    Q <- t(as.matrix(s1[, ..beta_cols]))
    SE <- t(as.matrix(s1[, ..se_cols]))
    snp_ids <- if ("snp_id" %in% names(s1)) s1$snp_id else paste0("snp", seq_len(ncol(Q)))

    return(list(taus = taus, Q_slope = Q, SE_tau = SE, snp_ids = snp_ids))
  }

  x <- readRDS(stage1_rds)
  if (is.null(x$taus) || is.null(x$Q_slope) || is.null(x$SE_tau)) {
    stop("stage1-rds must include taus, Q_slope, SE_tau")
  }
  Q <- as.matrix(x$Q_slope)
  SE <- as.matrix(x$SE_tau)
  if (nrow(Q) != length(x$taus) || nrow(SE) != length(x$taus) || any(dim(Q) != dim(SE))) {
    stop("Stage1 dimensions inconsistent: expected Q_slope and SE_tau as T x P")
  }
  snp_ids <- if (!is.null(x$snp_ids)) x$snp_ids else paste0("snp", seq_len(ncol(Q)))
  list(taus = as.numeric(x$taus), Q_slope = Q, SE_tau = SE, snp_ids = snp_ids)
}

parse_weights_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) stop("--weights is required")
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  out <- vector("list", length(parts))
  names(out) <- rep("", length(parts))
  for (i in seq_along(parts)) {
    kv <- strsplit(parts[i], "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Malformed weights entry: ", parts[i], " (expected name=path)")
    nm <- trimws(kv[1])
    fp <- trimws(kv[2])
    if (!nzchar(nm) || !nzchar(fp)) stop("Malformed weights entry: ", parts[i])
    out[[i]] <- fp
    names(out)[i] <- nm
  }
  out
}

read_weight_aligned <- function(path, taus_ref) {
  w <- readRDS(path)
  W <- if (is.list(w) && !is.null(w$W)) w$W else w
  w_taus <- if (is.list(w) && !is.null(w$taus)) as.numeric(w$taus) else NULL
  W <- as.matrix(W)

  if (!is.null(w_taus)) {
    if (nrow(W) != length(w_taus) && ncol(W) == length(w_taus)) W <- t(W)
    idx <- vapply(taus_ref, function(tau) {
      d <- abs(w_taus - tau)
      j <- which.min(d)
      if (d[j] > 1e-8) NA_integer_ else j
    }, integer(1))
    if (anyNA(idx)) stop("Weight taus do not align with stage1 taus: ", path)
    W <- W[idx, , drop = FALSE]
  } else {
    if (nrow(W) == length(taus_ref)) {
      # pass
    } else if (ncol(W) == length(taus_ref)) {
      W <- t(W)
    } else {
      stop("Weight matrix has no taus and dimensions do not match stage1 taus: ", path)
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
    ld <- determinant(diag(nrow(G)) + G, logarithm = TRUE)$modulus[1]
    val <- val + as.numeric(ld)
  }
  val
}

run_greedy_fixed_monotone <- function(taus, D_list, Kmax, method) {
  Kmax <- min(as.integer(Kmax), length(taus))
  selected <- integer(0)
  cur <- 0
  rows <- vector("list", Kmax)

  for (step in seq_len(Kmax)) {
    cand <- setdiff(seq_along(taus), selected)
    gains <- vapply(cand, function(ci) {
      newv <- obj_fixed_monotone(c(selected, ci), D_list)
      newv - cur
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

fit_wls <- function(y, se, W, idx_train, idx_hold) {
  y_tr <- y[idx_train]
  se_tr <- pmax(se[idx_train], 1e-10)
  X <- W[idx_train, , drop = FALSE] * (1 / se_tr)
  z <- y_tr / se_tr

  theta <- tryCatch(qr.solve(X, z), error = function(e) rep(NA_real_, ncol(W)))
  if (anyNA(theta)) {
    return(list(AIC = NA_real_, Q_train = NA_real_, Q_train_pt = NA_real_, Q_hold = NA_real_, Q_hold_pt = NA_real_))
  }

  pred_tr <- as.numeric(W[idx_train, , drop = FALSE] %*% theta)
  Q_train <- sum(((y_tr - pred_tr) / se_tr)^2)
  AIC <- Q_train + 2 * ncol(W)

  if (length(idx_hold) == 0L) {
    return(list(AIC = AIC, Q_train = Q_train, Q_train_pt = Q_train / length(idx_train), Q_hold = NA_real_, Q_hold_pt = NA_real_))
  }

  y_ho <- y[idx_hold]
  se_ho <- pmax(se[idx_hold], 1e-10)
  pred_ho <- as.numeric(W[idx_hold, , drop = FALSE] %*% theta)
  Q_hold <- sum(((y_ho - pred_ho) / se_ho)^2)

  list(
    AIC = AIC,
    Q_train = Q_train,
    Q_train_pt = Q_train / length(idx_train),
    Q_hold = Q_hold,
    Q_hold_pt = Q_hold / length(idx_hold)
  )
}

score_snp_models <- function(y, se, W_list, idx_train, idx_hold) {
  mods <- names(W_list)
  fits <- lapply(mods, function(m) fit_wls(y, se, W_list[[m]], idx_train, idx_hold))
  names(fits) <- mods
  aics <- vapply(fits, function(x) x$AIC, numeric(1))
  ord <- order(aics)

  best <- mods[ord[1]]
  list(
    best = best,
    delta = if (length(ord) >= 2) aics[ord[2]] - aics[ord[1]] else NA_real_,
    q_train_pt = fits[[best]]$Q_train_pt,
    q_hold_pt = fits[[best]]$Q_hold_pt,
    aics = aics
  )
}

run_model_compare_for_greedy <- function(greedy, taus, Q, SE, snp_ids, W_list, k_eval) {
  k_eval <- sort(unique(pmin(pmax(as.integer(k_eval), 1L), nrow(greedy))))

  full_idx <- seq_along(taus)
  full_best <- character(ncol(Q))
  for (i in seq_len(ncol(Q))) {
    sc <- score_snp_models(Q[, i], SE[, i], W_list, full_idx, integer(0))
    full_best[i] <- sc$best
  }

  detail_rows <- vector("list", length(k_eval) * ncol(Q))
  summary_rows <- vector("list", length(k_eval))
  winner_rows <- vector("list", length(k_eval))

  z <- 1L
  for (kk in seq_along(k_eval)) {
    k <- k_eval[kk]
    idx_train <- greedy$add_idx[seq_len(k)]
    idx_hold <- setdiff(seq_along(taus), idx_train)

    best_k <- character(ncol(Q))
    delta_k <- numeric(ncol(Q))
    qtrain_k <- numeric(ncol(Q))
    qhold_k <- numeric(ncol(Q))

    for (i in seq_len(ncol(Q))) {
      sc <- score_snp_models(Q[, i], SE[, i], W_list, idx_train, idx_hold)
      best_k[i] <- sc$best
      delta_k[i] <- sc$delta
      qtrain_k[i] <- sc$q_train_pt
      qhold_k[i] <- sc$q_hold_pt

      detail_rows[[z]] <- data.table(
        k = k,
        snp_id = snp_ids[i],
        winner = sc$best,
        winner_full = full_best[i],
        agree_full = as.integer(sc$best == full_best[i]),
        delta_train = sc$delta,
        q_train_pt = sc$q_train_pt,
        q_hold_pt = sc$q_hold_pt,
        overfit_gap = sc$q_hold_pt - sc$q_train_pt
      )
      z <- z + 1L
    }

    dt_k <- data.table(
      k = k,
      winner = best_k,
      agree_full = as.integer(best_k == full_best),
      delta_train = delta_k,
      q_train_pt = qtrain_k,
      q_hold_pt = qhold_k,
      overfit_gap = qhold_k - qtrain_k
    )

    summary_rows[[kk]] <- dt_k[, .(
      k = k[1],
      n_snps = .N,
      pct_agree_full = 100 * mean(agree_full, na.rm = TRUE),
      median_delta_train = median(delta_train, na.rm = TRUE),
      q25_delta_train = as.numeric(quantile(delta_train, 0.25, na.rm = TRUE)),
      q75_delta_train = as.numeric(quantile(delta_train, 0.75, na.rm = TRUE)),
      median_overfit_gap = median(overfit_gap, na.rm = TRUE),
      pct_overfit_gap_gt0 = 100 * mean(overfit_gap > 0, na.rm = TRUE)
    )]

    winner_rows[[kk]] <- dt_k[, .N, by = .(k, winner)]
  }

  list(
    detail = rbindlist(detail_rows, fill = TRUE),
    summary = rbindlist(summary_rows, fill = TRUE)[order(k)],
    winners = rbindlist(winner_rows, fill = TRUE)[order(k, -N)]
  )
}

save_compare_outputs <- function(out_dir, mode_name, greedy, res, greedy_all = NULL) {
  fwrite(greedy, file.path(out_dir, paste0("greedy_fixed_monotone_", mode_name, "_sequence.tsv")), sep = "\t")
  fwrite(res$summary, file.path(out_dir, paste0("model_compare_summary_by_k_", mode_name, ".tsv")), sep = "\t")
  fwrite(res$detail, file.path(out_dir, paste0("model_compare_detail_by_snp_", mode_name, ".tsv.gz")), sep = "\t")
  fwrite(res$winners, file.path(out_dir, paste0("winner_counts_by_k_", mode_name, ".tsv")), sep = "\t")

  p1 <- ggplot(greedy, aes(x = step, y = marginal_gain)) +
    geom_line(linewidth = 0.9, color = "#1b9e77") +
    geom_point(size = 1.3, color = "#1b9e77") +
    labs(
      title = paste0("Fixed-monotone greedy tau selection (", mode_name, ")"),
      x = "Step",
      y = "Marginal gain"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, paste0("greedy_marginal_gain_", mode_name, ".png")), p1, width = 8.5, height = 4.5, dpi = 170)

  p2 <- ggplot(greedy, aes(x = step, y = cum_frac)) +
    geom_line(linewidth = 0.9, color = "#d95f02") +
    geom_point(size = 1.3, color = "#d95f02") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste0("Fixed-monotone greedy tau selection (", mode_name, ")"),
      x = "Step",
      y = "Cumulative objective fraction"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, paste0("greedy_cumfrac_", mode_name, ".png")), p2, width = 8.5, height = 4.5, dpi = 170)

  p3 <- ggplot(res$summary, aes(x = k, y = pct_agree_full)) +
    geom_line(linewidth = 0.9, color = "#1f78b4") +
    geom_point(size = 1.5, color = "#1f78b4") +
    labs(
      title = paste0("Winner agreement with full tau grid (", mode_name, ")"),
      x = "Number of selected taus (k)",
      y = "Agreement with full grid (%)"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, paste0("agreement_vs_k_", mode_name, ".png")), p3, width = 8.5, height = 4.5, dpi = 170)

  p4 <- ggplot(res$winners, aes(x = factor(k), y = N, fill = winner)) +
    geom_col(position = "fill") +
    labs(
      title = paste0("Model winner composition by selected tau count (", mode_name, ")"),
      x = "k",
      y = "Proportion"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, paste0("winner_composition_by_k_", mode_name, ".png")), p4, width = 8.5, height = 4.5, dpi = 170)

  p5 <- ggplot(res$detail, aes(x = factor(k), y = delta_train, fill = factor(k))) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
    labs(
      title = paste0("Winner-vs-runner-up delta by selected tau count (", mode_name, ")"),
      x = "k",
      y = "Delta AIC"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, paste0("delta_boxplot_by_k_", mode_name, ".png")), p5, width = 8.5, height = 4.5, dpi = 170)

  if (!is.null(greedy_all)) {
    pcmp <- ggplot(greedy_all, aes(x = step, y = cum_frac, color = method)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.1) +
      labs(
        title = "Greedy cumulative objective: no_cov vs emp_cov",
        x = "Step",
        y = "Cumulative objective fraction",
        color = "Mode"
      ) +
      theme_minimal(base_size = 12)
    ggsave(file.path(out_dir, "greedy_cumfrac_compare_nocov_vs_empcov.png"), pcmp, width = 8.5, height = 4.5, dpi = 170)
  }
}

main <- function() {
  opt <- parse_args()
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$seed)

  if (!opt$greedy_mode %in% c("no_cov", "emp_cov", "both")) {
    stop("--greedy-mode must be one of: no_cov, emp_cov, both")
  }

  stage1 <- read_stage1(opt$stage1_tsv, opt$stage1_rds)
  taus <- stage1$taus
  Q <- stage1$Q_slope
  SE <- stage1$SE_tau
  snp_ids <- stage1$snp_ids

  if (!is.na(opt$max_snps) && opt$max_snps > 0 && ncol(Q) > opt$max_snps) {
    pick <- sample.int(ncol(Q), opt$max_snps)
    Q <- Q[, pick, drop = FALSE]
    SE <- SE[, pick, drop = FALSE]
    snp_ids <- snp_ids[pick]
  }

  w_map <- parse_weights_arg(opt$weights)
  W_list <- lapply(w_map, read_weight_aligned, taus_ref = taus)
  if (length(W_list) < 2) stop("Need at least two weight matrices for model comparison")

  k_eval <- opt$k_eval
  if (is.null(k_eval) || length(k_eval) == 0L) {
    k_eval <- sort(unique(c(9L, 13L, 17L, 25L, 33L, opt$kmax)))
  }

  greedy_results <- list()

  if (opt$greedy_mode %in% c("no_cov", "both")) {
    D_nocov <- build_contrast_rows(W_list)
    greedy_results$no_cov <- run_greedy_fixed_monotone(taus, D_nocov, opt$kmax, method = "no_cov")
  }

  if (opt$greedy_mode %in% c("emp_cov", "both")) {
    R_emp <- stats::cor(t(Q), use = "pairwise.complete.obs")
    eig <- eigen(R_emp, symmetric = TRUE)
    eig$values[eig$values < 1e-10] <- 1e-10
    R_emp <- eig$vectors %*% (eig$values * t(eig$vectors))
    D_empcov <- build_contrast_rows_empcov(W_list, R_emp)
    greedy_results$emp_cov <- run_greedy_fixed_monotone(taus, D_empcov, opt$kmax, method = "emp_cov")
  }

  if (opt$greedy_mode == "both") {
    cmp <- merge(
      greedy_results$no_cov[, .(step, no_cov_tau = add_tau)],
      greedy_results$emp_cov[, .(step, emp_cov_tau = add_tau)],
      by = "step", all = TRUE
    )
    fwrite(cmp, file.path(opt$out_dir, "greedy_nocov_vs_empcov.tsv"), sep = "\t")

    g_all <- rbindlist(greedy_results, use.names = TRUE)

    res_nocov <- run_model_compare_for_greedy(greedy_results$no_cov, taus, Q, SE, snp_ids, W_list, k_eval)
    res_empcov <- run_model_compare_for_greedy(greedy_results$emp_cov, taus, Q, SE, snp_ids, W_list, k_eval)

    save_compare_outputs(opt$out_dir, "no_cov", greedy_results$no_cov, res_nocov, greedy_all = g_all)
    save_compare_outputs(opt$out_dir, "emp_cov", greedy_results$emp_cov, res_empcov, greedy_all = g_all)

    message("Saved outputs in: ", opt$out_dir)
    message("no_cov summary:")
    print(res_nocov$summary)
    message("emp_cov summary:")
    print(res_empcov$summary)
  } else {
    mode_name <- opt$greedy_mode
    greedy <- greedy_results[[mode_name]]
    res <- run_model_compare_for_greedy(greedy, taus, Q, SE, snp_ids, W_list, k_eval)
    save_compare_outputs(opt$out_dir, mode_name, greedy, res)
    message("Saved outputs in: ", opt$out_dir)
    print(res$summary)
  }
}

main()
