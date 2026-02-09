#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

usage <- function() {
  cat(
"Usage:\n",
"  Rscript stage2/scripts/simulate_taugrid_model_recovery.R \\\n",
"    --weights add=path/add.rds,shift=path/shift.rds,both=path/both.rds,vqtl=path/vqtl.rds \\\n",
"    --out-dir results/sim_taugrid \\\n",
"    --n-snps 400 --n-reps 10 --k-eval 9,13,17,25,33,40 --grid-mode both\\n\\n",
"Optional args:\n",
"  --truth-probs add:0.25,shift:0.25,both:0.25,vqtl:0.25\n",
"  --effect-sd 0.08\n",
"  --noise-sd 0.04\n",
"  --rho 0.95\n",
"  --cov-jitter 0.25\n",
"  --seed 6520\n",
"  --help\n",
sep = "")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    weights = NULL,
    out_dir = "results/sim_taugrid",
    n_snps = 400L,
    n_reps = 10L,
    k_eval = c(9L, 13L, 17L, 25L, 33L, 40L),
    grid_mode = "both",
    truth_probs = "add:0.25,shift:0.25,both:0.25,vqtl:0.25",
    effect_sd = 0.08,
    noise_sd = 0.04,
    rho = 0.95,
    cov_jitter = 0.25,
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

    if (key == "weights") out$weights <- val
    else if (key == "out-dir") out$out_dir <- val
    else if (key == "n-snps") out$n_snps <- as.integer(val)
    else if (key == "n-reps") out$n_reps <- as.integer(val)
    else if (key == "k-eval") out$k_eval <- as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
    else if (key == "grid-mode") out$grid_mode <- trimws(tolower(val))
    else if (key == "truth-probs") out$truth_probs <- val
    else if (key == "effect-sd") out$effect_sd <- as.numeric(val)
    else if (key == "noise-sd") out$noise_sd <- as.numeric(val)
    else if (key == "rho") out$rho <- as.numeric(val)
    else if (key == "cov-jitter") out$cov_jitter <- as.numeric(val)
    else if (key == "seed") out$seed <- as.integer(val)
    else stop("Unknown argument: --", key)
  }

  out
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

read_weight_any <- function(path) {
  w <- readRDS(path)
  W <- if (is.list(w) && !is.null(w$W)) w$W else w
  taus <- if (is.list(w) && !is.null(w$taus)) as.numeric(w$taus) else NULL
  W <- as.matrix(W)
  if (!is.null(taus) && nrow(W) != length(taus) && ncol(W) == length(taus)) W <- t(W)
  list(W = W, taus = taus)
}

align_weight_list <- function(weight_map) {
  w_objs <- lapply(weight_map, read_weight_any)
  taus_ref <- w_objs[[1]]$taus
  if (is.null(taus_ref)) stop("Weights must include taus for this simulation script")

  out <- list()
  for (nm in names(w_objs)) {
    w <- w_objs[[nm]]
    if (is.null(w$taus)) stop("Weight file missing taus: ", nm)
    idx <- vapply(taus_ref, function(tau) {
      d <- abs(w$taus - tau)
      j <- which.min(d)
      if (d[j] > 1e-8) NA_integer_ else j
    }, integer(1))
    if (anyNA(idx)) stop("Tau mismatch while aligning weights for model: ", nm)
    out[[nm]] <- w$W[idx, , drop = FALSE]
    storage.mode(out[[nm]]) <- "double"
  }

  list(taus = taus_ref, W_list = out)
}

parse_truth_probs <- function(x, models) {
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  out <- setNames(rep(0, length(models)), models)
  for (p in parts) {
    kv <- strsplit(trimws(p), ":", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Malformed truth-probs entry: ", p)
    nm <- trimws(kv[1])
    pr <- as.numeric(trimws(kv[2]))
    if (!nm %in% models) stop("Unknown model in truth-probs: ", nm)
    if (!is.finite(pr) || pr < 0) stop("Invalid probability for ", nm)
    out[nm] <- pr
  }
  s <- sum(out)
  if (s <= 0) stop("truth-probs sum must be > 0")
  out / s
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

run_greedy <- function(taus, D_list, Kmax, method) {
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

make_base_sigma <- function(Tt, rho, noise_sd) {
  idx <- seq_len(Tt)
  R <- outer(idx, idx, function(i, j) rho^abs(i - j))
  (noise_sd^2) * R
}

simulate_one_rep <- function(rep_id, taus, W_list, truth_probs, n_snps, effect_sd, sigma_base, cov_jitter) {
  models <- names(W_list)
  Tt <- length(taus)

  true_model <- sample(models, size = n_snps, replace = TRUE, prob = truth_probs)
  beta_mat <- matrix(NA_real_, nrow = Tt, ncol = n_snps)
  se_mat <- matrix(NA_real_, nrow = Tt, ncol = n_snps)
  Sigma_list <- vector("list", n_snps)

  idx <- seq_len(Tt)
  for (j in seq_len(n_snps)) {
    m <- true_model[j]
    W <- W_list[[m]]
    k <- ncol(W)

    theta <- rnorm(k, mean = 0, sd = effect_sd)
    if (k > 1) {
      theta[-1] <- theta[-1] + rnorm(k - 1, 0, effect_sd / 2)
    }

    # Per-SNP covariance variation to mimic heterogeneous uncertainty.
    scale_j <- exp(rnorm(1, 0, cov_jitter))
    rho_j <- max(0.7, min(0.999, 0.95 + rnorm(1, 0, 0.03)))
    R_j <- outer(idx, idx, function(i, ii) rho_j^abs(i - ii))
    Sigma_j <- scale_j * ((sigma_base + diag(diag(sigma_base))) / 2)
    Sigma_j <- (diag(diag(Sigma_j))^0.5) %*% R_j %*% (diag(diag(Sigma_j))^0.5)

    z <- as.numeric(t(chol(Sigma_j)) %*% rnorm(Tt))
    beta_mat[, j] <- as.numeric(W %*% theta) + z
    se_mat[, j] <- sqrt(pmax(diag(Sigma_j), 1e-12))
    Sigma_list[[j]] <- Sigma_j
  }

  list(
    rep = rep_id,
    beta_mat = beta_mat,
    se_mat = se_mat,
    Sigma_list = Sigma_list,
    true_model = true_model
  )
}

fit_dwls <- function(y, se, W) {
  se <- pmax(se, 1e-12)
  X <- W * (1 / se)
  z <- y / se
  theta <- tryCatch(qr.solve(X, z), error = function(e) rep(NA_real_, ncol(W)))
  if (anyNA(theta)) return(list(Q = NA_real_, AIC = NA_real_))
  pred <- as.numeric(W %*% theta)
  Q <- sum(((y - pred) / se)^2)
  list(Q = Q, AIC = Q + 2 * ncol(W))
}

fit_gls <- function(y, Sigma, W) {
  C <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(C)) return(list(Q = NA_real_, AIC = NA_real_))

  X <- backsolve(C, W, transpose = TRUE)
  z <- backsolve(C, y, transpose = TRUE)
  theta <- tryCatch(qr.solve(X, z), error = function(e) rep(NA_real_, ncol(W)))
  if (anyNA(theta)) return(list(Q = NA_real_, AIC = NA_real_))
  r <- y - as.numeric(W %*% theta)
  q <- tryCatch(as.numeric(crossprod(r, solve(Sigma, r))), error = function(e) NA_real_)
  list(Q = q, AIC = q + 2 * ncol(W))
}

score_models <- function(y, se, Sigma, W_list, idx_tau, se_mode) {
  mods <- names(W_list)
  aics <- rep(NA_real_, length(mods))
  names(aics) <- mods

  y_sub <- y[idx_tau]
  se_sub <- se[idx_tau]
  Sigma_sub <- Sigma[idx_tau, idx_tau, drop = FALSE]

  for (m in mods) {
    W <- W_list[[m]][idx_tau, , drop = FALSE]
    fit <- if (se_mode == "dwls") fit_dwls(y_sub, se_sub, W) else fit_gls(y_sub, Sigma_sub, W)
    aics[m] <- fit$AIC
  }

  ord <- order(aics)
  list(
    winner = names(aics)[ord[1]],
    delta = if (length(ord) >= 2) aics[ord[2]] - aics[ord[1]] else NA_real_,
    aics = aics
  )
}

main <- function() {
  opt <- parse_args()
  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$seed)

  if (!opt$grid_mode %in% c("no_cov", "emp_cov", "both")) {
    stop("--grid-mode must be one of: no_cov, emp_cov, both")
  }

  w_map <- parse_weights_arg(opt$weights)
  aligned <- align_weight_list(w_map)
  taus <- aligned$taus
  W_list <- aligned$W_list
  models <- names(W_list)
  truth_probs <- parse_truth_probs(opt$truth_probs, models)

  sigma_base <- make_base_sigma(length(taus), opt$rho, opt$noise_sd)

  # Tau-grid selectors are built from model geometry; empirical-cov selector uses
  # an expected tau-correlation proxy (AR1 from --rho).
  D_nocov <- build_contrast_rows(W_list)
  greedy_nocov <- run_greedy(taus, D_nocov, max(opt$k_eval), method = "no_cov")

  grid_defs <- list()
  if (opt$grid_mode %in% c("no_cov", "both")) {
    grid_defs$no_cov <- greedy_nocov
  }
  if (opt$grid_mode %in% c("emp_cov", "both")) {
    R_proxy <- outer(seq_along(taus), seq_along(taus), function(i, j) opt$rho^abs(i - j))
    D_empcov <- build_contrast_rows_empcov(W_list, R_proxy)
    grid_defs$emp_cov <- run_greedy(taus, D_empcov, max(opt$k_eval), method = "emp_cov")
  }

  if (opt$grid_mode == "both") {
    cmp <- merge(
      grid_defs$no_cov[, .(step, no_cov_tau = add_tau)],
      grid_defs$emp_cov[, .(step, emp_cov_tau = add_tau)],
      by = "step", all = TRUE
    )
    fwrite(cmp, file.path(opt$out_dir, "sim_greedy_nocov_vs_empcov.tsv"), sep = "\t")
  }

  all_details <- list()
  row_idx <- 1L

  for (rep_id in seq_len(opt$n_reps)) {
    sim <- simulate_one_rep(
      rep_id = rep_id,
      taus = taus,
      W_list = W_list,
      truth_probs = truth_probs,
      n_snps = opt$n_snps,
      effect_sd = opt$effect_sd,
      sigma_base = sigma_base,
      cov_jitter = opt$cov_jitter
    )

    for (grid_name in names(grid_defs)) {
      for (k in sort(unique(pmin(pmax(opt$k_eval, 1L), nrow(grid_defs[[grid_name]]))))) {
        idx_tau <- grid_defs[[grid_name]]$add_idx[seq_len(k)]

        for (se_mode in c("dwls", "tau_cov")) {
          for (j in seq_len(opt$n_snps)) {
            sc <- score_models(
              y = sim$beta_mat[, j],
              se = sim$se_mat[, j],
              Sigma = sim$Sigma_list[[j]],
              W_list = W_list,
              idx_tau = idx_tau,
              se_mode = se_mode
            )

            all_details[[row_idx]] <- data.table(
              rep = rep_id,
              snp = j,
              grid = grid_name,
              k = k,
              se_mode = se_mode,
              true_model = sim$true_model[j],
              winner = sc$winner,
              correct = as.integer(sc$winner == sim$true_model[j]),
              delta = sc$delta
            )
            row_idx <- row_idx + 1L
          }
        }
      }
    }
  }

  detail <- rbindlist(all_details, fill = TRUE)
  summary <- detail[, .(
    n = .N,
    accuracy = mean(correct, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    q25_delta = as.numeric(quantile(delta, 0.25, na.rm = TRUE)),
    q75_delta = as.numeric(quantile(delta, 0.75, na.rm = TRUE))
  ), by = .(grid, k, se_mode)][order(grid, k, se_mode)]

  by_truth <- detail[, .(
    n = .N,
    accuracy = mean(correct, na.rm = TRUE)
  ), by = .(grid, k, se_mode, true_model)][order(grid, k, se_mode, true_model)]

  confusion <- detail[, .N, by = .(grid, k, se_mode, true_model, winner)][order(grid, k, se_mode, true_model, -N)]

  fwrite(rbindlist(grid_defs, use.names = TRUE), file.path(opt$out_dir, "sim_greedy_sequences.tsv"), sep = "\t")
  fwrite(summary, file.path(opt$out_dir, "sim_summary.tsv"), sep = "\t")
  fwrite(by_truth, file.path(opt$out_dir, "sim_accuracy_by_truth.tsv"), sep = "\t")
  fwrite(confusion, file.path(opt$out_dir, "sim_confusion_counts.tsv"), sep = "\t")
  fwrite(detail, file.path(opt$out_dir, "sim_detail.tsv.gz"), sep = "\t")

  p1 <- ggplot(summary, aes(x = k, y = accuracy, color = se_mode)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    facet_wrap(~grid) +
    labs(
      title = "Known-truth model recovery vs number of taus",
      x = "Number of selected taus (k)",
      y = "Recovery accuracy",
      color = "SE mode"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(opt$out_dir, "sim_accuracy_vs_k.png"), p1, width = 9, height = 5, dpi = 170)

  p2 <- ggplot(by_truth, aes(x = k, y = accuracy, color = se_mode)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.1) +
    facet_grid(true_model ~ grid) +
    labs(
      title = "Recovery accuracy by true model",
      x = "Number of selected taus (k)",
      y = "Recovery accuracy",
      color = "SE mode"
    ) +
    theme_minimal(base_size = 11)
  ggsave(file.path(opt$out_dir, "sim_accuracy_by_truth_vs_k.png"), p2, width = 10, height = 7, dpi = 170)

  p3 <- ggplot(summary, aes(x = k, y = median_delta, color = se_mode)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.3) +
    facet_wrap(~grid) +
    labs(
      title = "Winner-vs-runner-up margin vs number of taus",
      x = "Number of selected taus (k)",
      y = "Median delta AIC",
      color = "SE mode"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(opt$out_dir, "sim_delta_vs_k.png"), p3, width = 9, height = 5, dpi = 170)

  message("Saved outputs in: ", opt$out_dir)
  print(summary)
}

main()
